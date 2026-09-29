# Changelog

Notable changes to vdjtools v2. Releases before 3.0.0 are recorded in the git tags
(`v2.5.0` … `v2.9.0`) and their commit history.

## 4.6.1 — 2026-09-29

Everything in 4.6.0, which **built every wheel and published none**: it carried a
`test_changelog_covers_every_release` that asserted the checkout has release tags.
`actions/checkout` clones without them, so the assertion fired inside a fixture and errored all
four Python jobs on both platforms — and `publish.yml` runs that suite as a precondition for
upload, so nothing shipped.

The test is **deleted**, not repaired. It was added unprompted during a changelog cleanup and was
never part of the request. The `fetch-depth: 0` added to both workflows to feed it is reverted, so
they are byte-identical to before.

The changelog *content* fixes from that cleanup stay, because those were the request: the missing
3.9.2 and 3.17.0 sections, the stray second `## Unreleased` folded into the 3.14.0 that shipped it,
3.0.0's date, and the NOT PUBLISHED marker on 4.4.0.

**4.6.0 is not on PyPI and its GitHub Release is withdrawn; the tag stays for the history.**

## 4.6.0 — 2026-09-29

### Added — `io.strip_allele_values`: the allele strip, resolved once per distinct call

`strip_allele` is an expression, so it runs its split / strip / regex / unique / sort / join once
per **row**. A repertoire has 10^5-10^6 rows and 10^1-10^2 distinct segment calls, so it was
answering the same handful of questions hundreds of thousands of times. Measured on `TRBV` calls
with 261 distinct values:

| rows | distinct | expression | per distinct | gain |
|---|---|---|---|---|
| 10,000 | 261 | 4.28 ms | 1.10 ms | 3.9x |
| 50,000 | 261 | 20.63 ms | 1.74 ms | 11.9x |
| 200,000 | 261 | 82.45 ms | 3.29 ms | **25.1x** |
| 1,000,000 | 261 | 431.53 ms | 13.32 ms | **32.4x** |

`strip_allele_values(series)` is built **from** the expression rather than beside it, so the two
cannot state the rule differently; only its domain changes. Eight call sites that already held a
column now use it — four in `signature.features` (which run per sample per locus) and four in
`preprocess.batch`.

**Two faster expressions were measured and rejected.** Replacing `list.eval` with whole-string
regexes is **100.4 ms against 91.5**, i.e. slower: the cost is the list machinery (split / unique /
sort / join), not the per-row sub-expression. Casting to `Categorical` first is 92.5 ms, because
polars decodes a categorical back to strings for these operations instead of working on its
dictionary. There is no faster expression — the win is doing the work 261 times instead of a
million, which is also why this is **not** a C++ candidate: a native pass over every row would
still be a pass over every row.

`strip_allele` is unchanged and stays the right call for a LazyFrame, a `group_by` aggregation or
a filter predicate. The twelve call sites of that kind were left alone.

### Added — `io.translate_junctions`, the batched translation that existed and could not be reached

`_core.translate_junctions` shipped in 4.5.0 and made the format readers 5.32x, but the only
public translation was the scalar `io.convert.translate`: the C++ was reachable through
`vdjtools._core` or a private helper and nowhere else. A caller holding a nucleotide column had
no batched entry point. `vdjtools.io.translate_junctions(nt, threads=0)` is that entry point
— a `pl.Series` or any iterable in, a `pl.Series` out, input order, nulls preserved.

**66.5 -> 5.0 ms on 50,000 junctions of 39 nt, 13.4x**, and identical to
`to_unified_cdr3aa(translate(nt))` on every frame offset at every length from 3 to 59, with and
without non-ACGT bases.

**The scalar `translate` is not deprecated and is not the slow path for its own callers.** The
anchor checks in `model.collapse` and `model.reference` translate a **single codon**, and a
pybind11 round trip costs **498 ns against the dict walk's 264 ns** — routing those through C++
would be 1.9x slower. Batch when you hold a column; keep the dict when you hold a codon. That
measurement is why they were left alone rather than swept up.

`seqtree` has no translation and deliberately gets none: it carries `Alphabet::{AminoAcid,
Nucleotide, NucleotideIUPAC}` for indexing and scoring and nothing that maps codons to residues.
The genetic code belongs where the AIRR data model is, and two owners for it is how they drift.

### Python 3.10 had no test coverage at all for five CI runs, and the red X was read as a flake

`test_working_notes_stay_private.py` imported `tomllib`, which is stdlib only from **3.11**, while
`requires-python` is **>=3.10**. Because it fails at *import*, pytest reported
`Interrupted: 1 error during collection` and ran **nothing** on 3.10 -- `19 deselected, 1 error` --
so from `19ecbb8` (2026-09-28 15:54 UTC) the 3.10 leg of CI on both ubuntu and macos was a red X
with no signal behind it rather than one failing assertion. A collection error is not one broken
test; it is the whole suite quietly not running, and it looks identical to a flake in the job list.

Fixed with a conditional import and an explicit `tomli>=1.1; python_version<'3.11'` in the `[test]`
extra. `tomli` is already present wherever pytest runs on 3.10 -- pytest requires it there -- but
leaning on another project's transitive pin to keep a CI leg alive is the same shape as the bug.

**Runtime 3.10 support was never affected.** Verified from the *published* 4.5.0 `cp310` wheel in a
clean 3.10.20 venv rather than from this tree: `infer_nt`, `infer_nt_batch` (a declined row coming
back null), `generate(engine="native")`, `native.pgen_nt_batch` and `_core.translate_junctions` all
behave. The full suite on 3.10 is **1175 passed, 20 skipped**; the nine fewer passes than 3.12's
1184 are scirpy/dandelion skips in that environment, not a version defect.

`tomllib` was the only 3.11+ dependency in the tree -- `StrEnum`, `typing.Self`, `TaskGroup`,
`except*`, `datetime.UTC` and `itertools.batched` appear nowhere in `python/` or `tests/`.

## 4.5.0 — 2026-09-29

### Every per-item loop in the package, audited and fixed where it mattered

All **84 modules** were read for one pattern: Python iterating over data where a batched native
call, a polars expression, or a grouped computation does the same work. Twenty-one have no
Python-level iteration at all; most of the rest iterate over *structure* — columns, loci, genes,
events, marginals — which is not the same thing. What was found, and what was done:

**Two native entry points that did not exist.**

`native.pgen_nt_batch` is the batch nucleotide Pgen. It was the one native Pgen with no batch, so
any caller with more than one sequence either looped in Python or wrapped that loop in a
`ThreadPoolExecutor` handing out **one task per sequence** — the only pool in the package
dispatching per item rather than per contiguous slice, and the shape rejected everywhere else here.
Three call sites paid for it and all three are one call now: `model.score._pgen_nt_many` (its pool
is gone), the `pgen` CLI command (which splits nucleotide from amino-acid input once and makes one
batched call per kind), and `sc.paired_pgen`, which scores each chain's **distinct** clonotype keys
through `pgen_aa_batch` instead of looping `pgen_aa` per row.

`_core.translate_junctions` translates a whole junction column the way the legacy converters do —
including the bidirectional out-of-frame walk polars cannot express — and is checked against this
repo's own `to_unified_cdr3aa(translate(nt))` on 4,009 sequences at every frame offset, with and
without non-ACGT bases: **0 mismatches**.

**`generate(engine="native")` — the ancestral sampler in C++.** Measured on one core:

| locus | reference seq/s | native seq/s | gain |
|---|---|---|---|
| TRB | 25,428 | 6,196,346 | **244x** |
| TRA | 46,125 | 7,263,921 | 158x |
| TRG | 52,221 | 7,794,484 | 149x |
| IGH, `productive_only` | 4,029 | 2,325,649 | **577x** |
| TRD, `productive_only` | 6,057 | 3,373,886 | 557x |

IGH is the one that mattered: pool generation, not featurisation, is the dominant serial cost of a
synthetic corpus build, at 77 s of single-core work per 5.2M sequences.

Every draw is seeded from `(seed, row)`, so the output depends on the seed and the row index alone
and **never on the thread count** — identical at 1, 8 and auto. A corpus whose contents depend on
the builder's core count cannot be compared with one built anywhere else.

**The default engine stays `"reference"`.** The two are different random streams: every shipped
artifact drawn from this function was built with the Python sampler, so its `seed` stream is
frozen. They agree in *distribution* — against a 200,000-draw native sample, TRB V usage
correlation **0.99020 / 0.99406 / 0.99766** and total variation **0.0337 / 0.0247 / 0.0156** for a
reference draw of 4,000 / 10,000 / 20,000. It tightens as `1/sqrt(n_ref)`, which is what two
samplers of the same distribution do and what a mis-indexed deletion would not; nucleotide length
means agree to 0.036 nt and the productive fraction to 0.002.

**The format readers stopped looping over rows.** Nine converters each built one Python dict per
input row and handed `_finalize` a list of them, so a 42,877-row immunoSEQ export was 42,877 dict
allocations plus a schema inference pass over them before any of the work polars exists for. They
are expressions now, and **all ten shipped fixtures are byte-identical to the previous output**,
compared as written parquet — which is how the rewrite was driven, baseline first.

| fixture | rows | before | after | gain |
|---|---|---|---|---|
| immunoseq | 42,877 | 356.1 ms | 66.9 ms | **5.32x** |
| imgthighvquest | 7,199 | 117.8 ms | 31.5 ms | 3.74x |
| migec | 2,420 | 7.7 ms | 4.9 ms | 1.58x |
| rtcr | 693 | 5.5 ms | 4.3 ms | 1.29x |
| mixcr | 262 | 3.5 ms | 3.5 ms | 1.01x |

The small files are fixed overhead, not row cost — 8.31 to 1.56 microseconds per row is the change,
and it shows up where there are rows.

**The single-cell per-cell loops are window functions.** `sc.resolve_chains` and `sc.pair_chains`
each looped over cells, slicing a frame and calling `.to_dicts()` per cell — so the cost scaled
with the cell count rather than with the data. Output equality is asserted at every size:

| cells | contigs | `resolve_chains` | `pair_chains` |
|---|---|---|---|
| 400 | 929 | 267.3 → 5.1 ms (52x) | 495.3 → 8.8 ms (57x) |
| 4,000 | 9,536 | 2,694.4 → 8.7 ms (311x) | 4,999.4 → 18.2 ms (274x) |
| 20,000 | 47,657 | 12,690.5 → 24.6 ms (**517x**) | 24,924.7 → 46.1 ms (**540x**) |

The loop's row order is reproduced exactly — cells in first-appearance order, and within a cell
heavy then light then b_light, each by rank. Sorting by cell id instead would have silently
reordered every caller's rows.

**The second pass found the one that mattered most: `_align_init`.** EM seeds gene usage from an
alignment vote before the first iteration — each read voting its longest-matching V and J germline
— and that was a Python loop doing `nV + nJ` germline comparisons per read, each a character walk.
It runs before **every** fit, the native EM path included, and it also built a whole prepared model
it did not need. Now one native pass, with the votes still accumulated in read order so the seeded
tables are **bitwise-identical**:

| locus | genes | 20,000 reads, was | now | gain |
|---|---|---|---|---|
| TRB | 59 V x 13 J | 436.9 ms | 13.98 ms | **31.2x** |
| TRA | 47 V x 61 J | 572.1 ms | 15.61 ms | **36.6x** |
| IGH | 75 V x 6 J | 403.6 ms | 15.39 ms | 26.2x |
| TRG | 9 V x 5 J | 122.2 ms | 8.46 ms | 14.4x |

Splitting a tied vote is load-bearing and survives unchanged: germline-identical paralogs
(TRBV6-2/6-5/6-6, IGKV2-28/2D-28) tie exactly, and handing the family to one representative seeds
the rest at `P(V) = 0`, which the E-step's zero-probability skip then makes absorbing.

`biomarker.cooccurrence._dense` mapped sample ids through a Python dict per row of a
(features x samples) frame; that is one `replace_strict` pass. `collapse._rep` was a
`map_elements` at eight sites — **54,972 Python calls per IGH model load** for a string
concatenation — and is now an expression; the output is identical across all seven loci, and the
**wall time did not move** (0.90x-1.09x), because a model load is dominated by polars' own query
overhead rather than by those calls. It is kept as the better code, not as a speedup.

**Two more, measured.** `dynamics.expansion._exact_p` evaluated `betabinom.pmf` once per clonotype
when the pmf depends only on the total, and a repertoire's totals repeat heavily — 20,000
clonotypes over totals 1..60 carry 59 distinct ones. Grouped by total, with the
minimum-likelihood region found by position in the sorted pmf rather than a mask per row:
**10.3x to 307.5x** (614.6 → 3.31 ms at n=20,000). The region is identical; its mass is now summed
in ascending-probability order, so a p-value can differ in its last bits — max relative difference
**1.08e-14** over four cohort shapes and three dispersions, and ascending order is the more
accurate of the two. `model.data.write_prepared` builds the arda hand-off as one column instead of
a Python string per clonotype: **4.6x** (20.9 → 4.5 ms on 20,000 rows), byte-identical.

**Measured and deliberately left alone, with the numbers, so nobody re-derives this.** Steady
state on a 16-core M-series, after warm-up — which matters, because several of these look an order
of magnitude worse on a first call and that is polars query compilation, not the code:

| thing | measured | why it stays |
|---|---|---|
| `load_bundled` | 36 ms TRB, 57 ms IGH | parquet read plus `repair_anchors`; caching it is not allowed |
| `collapse_alleles` | 13-53 ms per locus | once per model load, over a few hundred alleles |
| `check_model` | 12.9 ms TRB, 24.4 ms IGH | a diagnostic, run on demand |
| `from_arda` | 9.3 ms | its two `map_elements` run over ~500 germline alleles, once |
| `prepare` | 10-15 ms cold | the reference path's tables |
| `entropy_table` / `mutual_information` | 8.7-19.1 ms | see below |
| `stitch_frame` | 1.01 microseconds/row | its germline lookups were hoisted long ago |
| `efron_thisted` | ~210 operations | bounded by `max_depth=20`, independent of clonotype count |

`model.analyze`'s group-wise entropies were **converted and then reverted**: a first-call profile
said `entropy_table` was 386 ms and 90% of it was polars `group_by` iterated in Python, so the loops
were rewritten as aggregations. Bit-identical output, and **slower** — 0.6x on TRB
`mutual_information`. The steady state was 8.7 ms all along and aggregation overhead on a
few-hundred-row table exceeds what slicing costs. Recorded because the first measurement was the
misleading one, not the code.

`io.read_vidjil` keeps its loop because it parses a JSON document, not a table — there is no column
to express the work over — but it now declares its schema, so an all-null D call cannot come out as
dtype `Null`. `overlap.similarity`'s identity block and `biomarker.metaclonotype`'s union-find do
work proportional to their output, not to the input, and `overlap.alice` already makes one
`pgen_aa_batch` call per allele pair.

**Every pool in the package was re-checked** and each wraps something that releases the GIL or is a
process pool over contiguous slices: `io.batch.map_samples` (polars read, 3.85x at 8 workers),
`model.data.build_all` (arda subprocess + native E-step), `signature.cohort.parallel_rows` and
`signature.corpus`'s two builders (spawned processes, no serial fallback). The one that dispatched
per item is gone.

### The thread-scaling guard needed two bars, not one

`test_quadrupling_the_threads_roughly_quarters_the_wall_time` demanded 2.5x at 4 threads; a 4-vCPU
CI runner measures **2.25x**, which is real parallelism at 56% efficiency and not a defect. It is
two bars now: **1.8x at 4 threads everywhere** (a GIL-bound stage gives ~1.7x there, and Amdahl on
a 35% serial fraction caps it at 1.96x even with no contention), and **2.5x only where 8 or more
cores are usable**, which is where near-linear scaling is visible and where the 1.86x the Python
codon reconstruction used to give would be caught. `available_cores` is the gate, not
`os.cpu_count` — a CI container's quota is invisible to the latter.

This is why **4.4.0 built every wheel and published none**: the guard failed the pre-publish test
job. The library code in it is what ships here.

## 4.4.0 — 2026-09-29

**NOT PUBLISHED.** Its own pre-publish thread-scaling guard failed on a 4-vCPU runner, so
every wheel built and none was uploaded; the tag stays for the history and the GitHub
Release was withdrawn. Everything below shipped in **4.5.0**, whose section records the
corrected guard. Nothing needs to be installed from here.

### `infer_nt` batched, and the codon reconstruction moved into C++ (#181)

`infer_nt` wrapped a native DP in per-row Python, and on a real annotation table that wrapper was
the cost, not the DP: **156.19 s of a 181.89 s stage, 85.9%**, at ~200,000 distinct junctions in
vdjdb-db's build. Two changes, in that order, because the first one measured why the second was
needed.

**`infer_nt_batch(model, cdr3_aas, v=, j=, n_best=8, threads=0)`** returns a frame with **one row
per input row in input order**, carrying `Scenario`'s fields as columns; a row the model cannot
explain is present with nulls, because the caller joins positionally. (Deliberately unlike
`native.best_aa_scenarios_batch`, which returns k rows per query and so signals a declined query by
absence.) `infer_nt` delegates to the same core, so the two cannot drift.

**Batching alone bought 1.34x, and the decomposition said why.** Stage 1 and the marginal re-score
were already native; laying a scenario out and picking its codons was Python, per candidate — 35%
of a human TRB row and **87% of a TRA** one. It held the GIL, so the batch was Amdahl-capped
whatever `threads` it was handed, and the TRA batch was flat in `threads` at 1.55x from 1 through
14. That is now `infer_nt_batch` in `src/pgen.cpp`: the same 25-state max-product DP, factor for
factor, including its tie-break — two paths of exactly equal weight resolve to the lexicographically
larger nucleotide string, so a run stays reproducible. The states carry no path (a backpointer table
does, and a position's nucleotide is read off its state), so the hot loop copies no strings.

With the search, the reconstruction and the re-score all in one call, the GIL is released for the
whole batch instead of reacquired once per row — which is what makes `threads` worth anything.

**Measured**, released VDJdb 2026-06-03 rebuilt into the distinct `(cdr3, v.segm, j.segm)` key set
the way the consumer builds it, samples of 3,000 from 111,655 allele-resolvable human TRB keys and
52,191 TRA, 16-core M-series (`threads=0` = 14 workers):

| | TRB ms/key | TRA ms/key |
|---|---|---|
| was: per-row `infer_nt`, serial | 1.7010 | 0.7938 |
| was: per-row `infer_nt` in a 4-thread slice pool | 0.8529 | 0.6682 |
| was: `infer_nt_batch`, Python reconstruction, `threads=0` | 0.6376 | 0.5080 |
| now: per-row `infer_nt`, serial | 1.0060 | 0.0791 |
| now: `infer_nt_batch`, `threads=0` | **0.0863** | **0.0071** |

**9.89x on TRB and 94.05x on TRA** against the threaded per-row loop; 10.30x / 98.58x at
`threads=16`. Over the whole human key set that is **130.1 s to 10.0 s, 13.01x**. Single-row
`infer_nt` is faster too — 1.69x on TRB and 10.04x on TRA — because the reconstruction it does once
per candidate is the same DP.

**Threads now scale, which is the whole test of whether a pool is parallelism or overhead with
extra steps.** Doubling the workers roughly halves the wall time: on TRB, 1.96x / 3.75x / 7.37x /
11.57x at 2 / 4 / 8 / 16 threads against 1. `tests/python/test_infer_nt_batch.py` asserts 2.5x at 4
threads, raised from the 1.25x that was honest while the reconstruction held the GIL. A caller must
**not** wrap this in a pool of its own: one batched call using the library's own kernel threads is
the supported shape, and a process pool around it oversubscribes the machine.

**It is the per-row loop's answer, not a new one** — every field of every row, at any thread count,
across all four V/J call modes, checked against the previous Python reconstruction on 4,134
row-comparisons over three loci: **0 mismatches**. The invariant test earned its place twice: the
first batch reported `n_candidates` *after* the top-`n_best` truncation, so a multi-allele call that
pooled 12 candidates reported 8 while every other field stayed correct; and the first native version
folded an **empty** candidate list into `None`, answering a call that named no allele by
marginalizing over all of them.

**A gene-level V/J call reached an answer for the first time.** The scenario search has always
resolved `TRBV18` to a representative allele; the marginal re-score did not — it was handed the
caller's raw name, and `pgen_nt` raises on a gene name on purpose (the `_gene_idx` trap), so a call
stage 1 accepted, stage 2 rejected with a `KeyError`. It hit 19 of 400 sampled human TRA keys. Every
name is now resolved once to a model index before either stage sees it.

**The `n_best` default stays at 8, and now there is less reason than ever to lower it.** It used to
be the knob that mattered because it multiplied the Python reconstruction; with the reconstruction
native it buys 1.39x on TRB and 1.19x on TRA, against a real cost in agreement:

| `n_best` | TRB ms/key | same `cdr3_nt` | TRB rows with no runner-up | TRA ms/key | same `cdr3_nt` | TRA rows with no runner-up |
|---|---|---|---|---|---|---|
| 8 | 0.0910 | 100% | 122 of 2,961 | 0.0076 | 100% | 135 of 2,711 |
| 4 | 0.0655 | 97.7% | 347 of 2,961 | 0.0064 | 99.9% | 433 of 2,711 |
| 2 | 0.0502 | 90.4% | 1,066 of 2,961 | 0.0059 | 99.2% | 1,359 of 2,711 |
| 1 | 0.0437 | 77.0% | 2,961 of 2,961 | 0.0055 | 92.9% | 2,711 of 2,711 |

Agreement is against the `n_best=8` default on the same 3,000-key samples. On the rows that do move
the alternative is not implausible — its marginal Pgen is a median 0.892 of the default's at
`n_best=4` — but `Scenario.margin` degrades well before `cdr3_nt` does, because a low `n_best`
leaves rows with no second candidate to compare against. `n_best=1` is the joint argmax with no
marginal re-score at all, and it has no `margin` by construction.

### Documentation

The documentation was surveyed page by page against the code it describes. What it got wrong:

- **Six of the fifteen example notebooks did not run.** Three referenced a column `cdr3_aa` that
  the canonical schema does not have — it is `junction_aa`, and the distinction is two residues of
  anchor, not a spelling. `aging.py` used `plt` without importing matplotlib anywhere and
  `emerson_cmv_hla.py` used `time` the same way; `model_explorer.py` bound `fig` in two cells,
  which marimo refuses outright. All fifteen now run end to end as plain scripts, and
  `tests/python/test_examples_wiring.py` checks the cell graph statically so this class cannot
  return.
- **`docs/cli.rst` documented 6 of the 10 signature options that change an answer**, omitting
  `--named` — the headline feature of 4.3.0 — along with `--winsor-p` and `--on-duplicate`, and
  33 model-workshop options had no prose anywhere. All are documented now, and
  `test_cli_help.py` fails if a flag that changes an answer is added without prose.
- **The corpus kinds were miscounted**: four are synthetic and five are real, not five and three.
- **`examples/README.md` told readers to run `vdjtools presets`** and `--preset classify`, both
  deleted in 4.0.0.
- `docs/notebooks.rst` listed 7 of the 15 notebooks; it now lists all of them, grouped by the
  question each one answers.

What it was missing, for a reader who is an immunologist rather than a statistician:

- Fourteen glossary entries for the statistical vocabulary the output uses — `channel`,
  `principal component`, `clr`, `logit`, `arcsine transform`, `robust z-score`, `support`,
  `named block`, `n_eff`, `dispersion`, `effective dimension`, `Rao quadratic entropy`,
  `clonality`, `isotype`, `SHM`, `Zipf law` — each written as what it is for, not as a definition.
- A plain-language account of the three signature stages ahead of the diagram that names them.
- An explicit warning that `pgen_aa` and `infer_nt` take a **junction** despite naming their
  argument `cdr3_aa`. Passing a trimmed CDR3 does not raise; it scores a different rearrangement.
- Worked-examples and glossary cards on both landing pages.

## 4.3.0 — 2026-09-28

Three things a caller needs, and none of them is a new statistic. Everything below already existed
inside the library; what changed is what you can reach. `ISSUES.md` item 13.

### Added: `named=` — the named statistics, from the call that computes them

The rotation is fitted on features that have domain names — diversity, depth, clone-size fractions,
junction length, isotype composition, SHM, cross-locus yield — and then emits coordinates that do
not. Those features were computed for every sample and thrown away, and the only route to a Hill
number was three calls, a hand-written block list, and knowing that `apply` drops its own input.

```python
vsig(sample, corpus, n_components=32, named=True)            # or named=("div", "depth")
vsig_cohort(samples, corpus, n_components=32, named=True)
```
```bash
vdjtools signature --corpus blood --components 32 --named all -o sig.tsv
```

`named=()` is the default and reproduces the previous output exactly. `mir.signature.rsig` takes the
same argument for `depth` / `band` / `band_igh`.

Which blocks are reportable is **declared**, on `RawGroup(named=True)`, not inferred from width:
`pchem` is 30 static columns and is not reportable, `shm` is one column and is. Values carry their
declared transform rather than a natural scale — a `log10` diversity comes back as `log10`, and a
`clr` composition has no unique inverse — so `vsig:div:TRB:1D_c = 2.386` is not a clone count.

**Why this is not merely convenient**: the diversity floor of any study using this library is made
of exactly these blocks. Without a documented route to them there is no floor, and a signature that
beats nothing gets reported as if it beat something.

### Added: `channel_table()` and a `kind` column — what can I get, and in what units

One call instead of three lookups and a reading of `corpus.py`. Every `(block, feature)` with its
transform, support, emitted loci, and a `kind`: `rotated` (you get PCs instead), `named` (also
available in its own right), `channel` (carried through untouched). `--describe` prints the same
three kinds, and now a `transform` column, for the exact columns an invocation will emit.

### Added: `fit_cohort` — fit a corpus on your own cohort

`corpus.fit` was public, correct, and documented only as the shipped-artifact route; it takes
feature **rows**, so a caller holding repertoires had to run the featuriser per sample and assemble
the dicts. `fit_cohort` takes frames, runs the existing process pool, and returns the same `Corpus`
type the published artifacts are — `save`/`load` round-trips, and a collaborator can score against
it exactly as against a shipped one.

```python
corpus = fit_cohort(train_samples, sig="vsig", name="my-cohort", n_components=64, n_jobs=0)
corpus.save("my-cohort.npz")
```

`sig="rsig"` resolves through a **featuriser registry** (`register_featuriser`) rather than an
import, because nothing in vdjtools may import `mir`.

**Documented with the caveat a caller most needs.** Measured on 5,376 raw columns over seven loci
with one logistic head, at matched widths: a rotation fitted on 612 training repertoires beat the
shipped `blood` artifact on the training out-of-fold number at **4 of 5** widths (median 0.5039
against 0.4968), and on an external cohort it led at **0 of 5** (median 0.4757 against 0.5270).
Fitting in-cohort improves the number the configuration is *selected* on and neither held-out
read-out — a rotation fitted on one trial learns that trial's covariance. Both routes are
legitimate; the docs now say which answers which question.

### Changed: `mask` is documented as a feature block, not as QC

`vsig:mask:<L>:present` / `:estimable` were listed beside `v_fallback_frac` under quality control,
so the natural reading was "diagnostics, drop before modelling". That reading costs real signal.
Measured on 874 labelled samples with a second cohort held out entirely, adding the five presence
flags to the recommended feature set moved the **external** ROC-AUC from **0.6243 to 0.6676** for 17
extra columns. In that cohort 720 of 7,399 (sample, locus) pairs sit below a five-clonotype floor,
almost all TRD and TRG, and which donors those are tracks lymphocyte content rather than noise.

`docs/channels.rst` now groups `mask` separately from `qc`, with the distinction stated: `qc`
describes how much to trust the row, `mask` describes the donor.

### Note: no corpus is invalidated

A channel is pass-through — the rotation is indexed by `raw_columns` and nothing else, and
`corpus.apply` fills every registered channel from the sample. The new `rsig` channel families
appear in the output of every artifact already published, with no refit and no new download.

## 4.2.0 — 2026-09-28

### Added: the three real corpora — `blood`, `tissue`, `deep-tcr`

The seven corpora the docs have named since 4.0.0 all exist now. The three real ones are fitted on
repertoires rather than drawn from the recombination models, which is the whole point of having both
kinds: a synthetic corpus spans a *measured* range by construction, a real one carries the joint
structure the generative model does not produce — selection-shaped V/J usage, isotype and SHM
structure, and the cross-locus covariance of libraries prepared together.

| corpus | samples | study_ids | fitted on |
|---|---:|---:|---|
| `blood` | 11,117 | 947 | public bulk RNA-seq blood, capped at 30 per study |
| `blood-uncapped` | 22,441 | 947 | the same population, uncapped |
| `tissue` | 21,131 | 1,934 | public bulk RNA-seq non-blood, capped at 30 per study |
| `tissue-uncapped` | 33,874 | 1,934 | the same population, uncapped |
| `deep-tcr` | 3,936 | 7 | targeted/amplicon deep TCR, TRA + TRB |

**The population had to be re-derived, not inherited.** The previous reference slice selected on
`has_TRB | has_IGH`, and that field carries two incompatible definitions — `>= 100 reads` in the base
metadata table and "any row present" in roughly 24,600 appended rows — so a population built on it
admits thousands of ultra-shallow repertoires. The harmonized sample table has no per-locus read
column at all, so the floor is now applied against reads **measured off the store** in one streaming
pass. Blood comes out at 22,441 task-disjoint samples against the 23,234 recorded for the old slice;
that difference is the defect being removed.

**TRG and TRD get their own floor, and it is a measurement.** They cannot support 100 reads at scale:
at that floor tissue TRG has 367 samples in 78 study_ids, and at 30 it has 1,574 in 317; blood TRG
goes from 5,686 to 13,381. So the five main loci keep 100 and TRG/TRD use 30, recorded per artifact.
Each locus is fitted only on the rows that observed it, so this is a stratum rather than a filter on
the corpus.

**Both cap variants ship.** The per-study cap is the dominance control the predicate specifies — the
top 10 submissions hold 17.7% of SRA blood — and shipping the uncapped fit beside it makes the cap's
effect measurable rather than asserted. Emission is the expensive half and is done once; each extra
fit is seconds.

A real corpus's manifest carries no cohort label, no dataset name and no accession. What it records is
the population as *rule text*, the per-locus floors, aggregate counts, the measured read band per
locus, and the germline fingerprint the V/J columns are indexed by. A rotation is about features and
components; a cohort label is not one of its inputs.

### Changed: k = 256, and every component ships

All nine corpora are fitted at **256 components** per locus, up from 128, and the full rotation ships.
Truncating downward at apply time is exact, so a 256-component artifact serves any narrower request
while the reverse needs a refit.

**Widening is a strict superset, measured rather than assumed**: refitting `synthetic-blood` from 128
to 256 reproduces all 128 original components with `|cos| = 1.000000` on TRB, IGH and TRG, and
identical eigenvalues to `0.000e+00`. Existing numbers at 128 do not move.

It buys real variance on the statistics half — `naive` TRB goes 0.4236 at 128 to 0.6283 at 256, and
`tissue` TRB reaches 0.8541 — while on the geometry half 256 is well past diminishing returns: `rsig`
reaches 0.90 in **7 to 18** components and is effectively full rank by 256, so `--components 32` is
usually the right request there.

### Changed: corpus artifacts are release assets, fetched on first use

At k=256 a `vsig` artifact is ~10 MB and the nine corpora are ~110 MB across both halves — too much
for every `pip install` to carry in order to use one of them. So the wheels now ship a few KB of
**index** instead: `corpora.json`, naming every corpus with its size and SHA-256. The artifacts are
published as GitHub release assets under their own `corpora-*` tag, and the library fetches one on
first use into `$VDJTOOLS_CORPUS_DIR` (default `~/.cache/vdjtools/signature`), verifying it against
the shipped digest before it is visible under its cache name.

The resolution order is a local path, then the wheel, then the cache, then the release — so a caller
who fitted their own corpus and passes its path is never served a download of the same name. A
download that fails its digest is deleted rather than cached: a half-written rotation that loads is
far worse than one that is missing. `vdjtools corpus --fetch all` / `mir corpus --fetch all` pre-warms
the cache, which is the install-time step.

`publish.yml` now publishes **only for a `v*` tag** in both repos. Without that guard, creating the
`corpora-*` release would fire the publish workflow and try to upload whatever version `pyproject`
happened to carry — a release that is not a release of the package. Every job is gated, not just the
upload step: building twelve wheels for a data release is waste even when the upload is skipped.

The release tag is deliberately not the library version. An artifact changes far less often than the
code, and keying the download on the version would invalidate every cached corpus on a patch release
— the same mistake the germline fingerprint already replaced for the load gate.

### Fixed: four amplicon samples were being counted twice, with their loci split

In `deep-tcr`, four samples have their TRA and TRB files in *different* source mirrors. Grouping by
sample within each cohort turned each into two half-samples — losing their cross-locus pair columns
and double-counting them — and the fitter's `sample_id` join then fanned 3,940 rows into 3,948,
mapping raw rows to the wrong samples while the shapes still lined up. Enumeration is now keyed by
sample across all cohorts, and the fitter asserts uniqueness on both sides of that join rather than
trusting it.

## 4.1.0 — 2026-09-28

### Added: `synthetic-blood` and `synthetic-tissue`, two corpora that describe a real compartment

Four corpora now ship, all of them synthetic. The two new ones draw every repertoire as a
naive/memory **mixture** across three quantile ladders measured per locus on the compartment they are
named after, so a corpus spans the depth and clone-size range the samples it will be applied to
actually have.

**The gap they close is an order of magnitude.** `naive` and `memory` draw their depth over
`DEPTH_SPREAD`, 2.4x on TRD to 11.0x on IGH, which came from 1,168 deep blood samples. Measured
across the whole harmonized AIRR store — one streaming pass, every sample with at least 100 reads in
the locus, grouped by `sample_id` — the real axis is far wider: blood TRB clonotype richness spans
**43x**, 74 to 3,162 clonotypes between its 5th and 95th percentiles over **34,365 samples**, and
tissue IGH spans **259x**, 40 to 10,352 over 47,031. Bounds, centre and per-PC scaling are all
estimated from the draw and none of them extrapolates, so a corpus drawn across a tenth of the axis
cannot standardise the rest of it.

Three ladders per locus, five quantiles each, in `corpus.COHORT`:

- **richness**, the clonotype count;
- **reads per expanded clone**, `(reads - singletons) / (richness - singletons)`;
- the **singleton fraction** `f1` — in a reasonably deep library this is what stands in for the naive
  compartment, and it falls with donor age (Britanova et al., J Immunol 2014,
  [10.4049/jimmunol.1302064](https://doi.org/10.4049/jimmunol.1302064); 2016,
  [10.4049/jimmunol.1600005](https://doi.org/10.4049/jimmunol.1600005)). Blood TRB spans 0.215 to
  0.949 with a median of 0.759.

V and J usage need no ladder: they come out of the rearrangement model.

**The mixture is constructed, not sampled.** `round(f1 * N)` clonotypes get one read and the rest
share the remainder with Zipf rank-abundance frequencies on a floor of 2, so a sample hits its drawn
richness, read count and singleton fraction exactly. Contrast `memory`, where the multinomial's zeros
are dropped and both the richness and the singleton fraction come out of the sampling rather than
being asked for — which is why a `memory` corpus cannot be pointed at a measured cohort.

**Reads per *expanded* clone rather than per clonotype, because the obvious quantity is not
drawable.** A repertoire with `f1` singletons whose other clones all carry at least 2 reads has at
least `2 - f1` reads per clonotype. That is an identity, and a pair drawn from those two marginals
lands in the region it forbids about half the time on blood TRB; drawing the read count itself is
worse, at **21.1%** of 20,000 blood TRB draws. Reads per expanded clone is >= 2 whatever `f1` is, so
the trio has no forbidden region and nothing has to be clipped back out of one.

**Five quantiles, not two, and the resolution is what makes the corpus land.** The read count is a
product of all three drawn quantities, so it is the sharpest check on the draw. Against blood TRB's
measured median of 763 reads, a 40,000-draw simulation gives **762** at five quantiles and 846 at
three; blood IGH's 1,298 comes out 1,276 against 1,490.

**The three are drawn through a Gaussian copula at the cohort's measured rank correlations, for the
rotation rather than for any marginal.** A corpus's rotation *is* its covariance structure, and `f1`
against expansion size is -0.503 on blood TRB and -0.846 on blood IGK (n = 43,678) — fewer
singletons, bigger expansions. Drawing independently would hand the PCA a correlation the cohort does
not have, which no amount of correct marginals repairs.

The correlation matrix is factorised **in closed form rather than by an eigendecomposition**, and that
is a reproducibility requirement. `numpy.linalg.eigh` returns eigenvectors whose sign is a LAPACK
convention: negating a column leaves the covariance, every marginal and every correlation untouched
while changing the realised sample, so a corpus built where LAPACK signs a column differently would
differ byte-for-byte from the same build here while being statistically identical — and a corpus whose
values depend on the builder's linear-algebra library cannot be compared with one built anywhere else.
A 3x3 Cholesky factor is unique and is six arithmetic operations.

Acceptance, on 500 drawn samples per cohort against the cohort's own percentiles: the singleton
fraction lands within 0.03, and richness within 10%, on **all 14** (cohort, locus) pairs; the read
count is within 15% on 10 of 14. The exceptions are the four shallow tissue TR loci, where the read
count runs 20-52% high (tissue TRB 467 against 308) — a product of three heavy-tailed factors has a
median above the product of their medians unless the joint tails match exactly, and no
marginals-plus-copula draw does that. The drawn p95 of the read count is also *lower* than the
cohort's (blood TRB 3,755 against 5,274), so the fitted depth ceiling is roughly 29% tighter than the
cohort's own p95. The full table is in `docs/signature.rst`.

`naive` and `memory` are unchanged and remain as the pure-regime references: the shipped artifacts
still rebuild bit-for-bit, checked by building a `memory` corpus against this commit and against
4.0.1 and comparing the bytes.

### Added: `--depth-spread`, and `--size auto`

A corpus describes only the depths it was drawn across, so the range is now a flag on both `corpus`
commands rather than a constant. `--size 3162 --depth-spread 1000` draws log-uniformly from 100 to
100,000 receptors per locus. On a `synthetic-*` corpus the band is the cohort's own measured one and
`--depth-spread` is **refused** rather than silently ignored; `--size` there rescales the whole ladder
about its median, so `--smoke` is the same cohort shape at a shallower nominal depth.

`--size` defaults to `auto`: 10,000 for a pure regime, the cohort's median richness per locus for a
`synthetic-*` one.

### Changed: `synthesize` takes a corpus name

The first parameter of `vdjtools.signature.corpus.synthesize` is `corpus_name`, not `regime`, and it
takes one of `SYNTHETIC` — `naive`, `memory`, `synthetic-blood`, `synthetic-tissue`. `corpus_plan`
resolves a name to `(regime, cohort, size, depth_spread)` and is the one place that mapping exists,
so the `vsig` and `rsig` halves cannot resolve the same name differently. `mixed` is the regime and
is refused as a corpus name, with a message saying why: a mixture has no unparameterised form.

Also new or changed, all in `signature.corpus`: `COHORT`, `COHORT_QS`, `SYNTHETIC`, `cohort_bands`,
`ladder_draw`, `copula_uniforms`, `draw_plan`, `corpus_meta` (one manifest builder for both halves);
`draw_sample(pool, size, regime, rng, frac=None, mexp=None)`; `draw_one(pools, plan, ...)` over the
per-sample plan rather than a sizes dict; `pool_target(..., cohort=None)`; `depth_spread_of` and
`resolved_size` accept a cohort name. `SINGLETON_FRAC`, a placeholder from an unreleased commit,
never shipped.

The manifest gains `cohort`, `richness_band`, `expanded_count_band`, `singleton_frac_band`,
`rank_corr`, `depth_spread` and `depth_spread_requested`, and the `rsig` half records the depth fields
for the first time — it previously carried none, which left the one thing a reader needs in order to
know whether a corpus covers their samples readable only from the `vsig` half.

## 4.0.1 — 2026-09-27

Everything in 4.0.0 below, plus the two things that kept it off PyPI. **4.0.0 never published** — the
tag was pushed but no GitHub Release was created, and `publish.yml` fires on `release: published`, not
on a tag. The `v4.0.0` tag stays where it is as the record; 4.0.1 is the first 4.x release to reach
PyPI.

### Fixed: apply-time truncation is exact on OpenBLAS too, not only on Accelerate

`test_apply_time_truncation_is_bitwise_exact` passed on macOS and failed on both Linux jobs, at
`vsig:pc:TRB:PC01` = -0.4881794592306904 against -0.4881794592306884 — a relative 4e-15, which is
BLAS reassociation rather than a different rotation. `z @ rotation[:, :k]` is the same arithmetic as
the untruncated `z @ rotation` in a *different* GEMM shape, and OpenBLAS is entitled to accumulate a
2-column product in a different order than a 128-column one.

The claim the test pins — fit once at `0.95`, apply later at a fixed count, and get bit-identical
values — is the whole reason apply-time truncation exists instead of a refit, so the fix is to make it
true by construction rather than to relax the test: rotate through the whole stored rotation and slice
the *scores*. Both paths now evaluate one identical expression. Cost is one full gemv per locus
instead of a truncated one, microseconds on a (2100, 128) rotation, once per sample.

### Changed: `arda-mapper>=2.30.1`

arda 2.30.0 fixed `cdr3fix` placing the amino-acid V/J boundary up to two residues too far into the
junction: `_align` broke a tie by consuming more query, so residues that scored nothing were credited
to the segment. Every `v_end` / `j_start` vdjtools reads out of an amino-acid annotation moves with
it. The floor is 2.30.1 rather than 2.30.0 because 2.30.0 never reached PyPI either.

## 4.0.0 — 2026-09-27

**Signatures rewritten from scratch.** No legacy path, no backward compatibility, no artifact
carried over. The statistics half only; the geometry half follows in mirpy 4.0.0.

### The defect this fixes

The old system fitted a rotation on **10,000 individual clonotypes** from a prototype panel, while
every one of the 399 PC columns that rotation produced was a **repertoire** statistic. Its centre
and scale were fitted on **zero rows** of that artifact and arrived from a separate corpus of real
samples. Two independent fits, stitched — which is how one shipped reference came to pair a centre of
exactly `0.0` with a scale plainly fitted from data, putting a corpus-typical sample **81 robust
deviations** out, and how `tissue` came to carry `blood`'s coverage constants for all seven loci to
17 significant digits while its own location and scale had genuinely been refitted.

Now everything comes out of **one pass over one matrix of repertoires**: bounds, centre, scale,
rotation and per-PC scaling, per locus, from the same corpus.

### Fixed, found while vectorizing the feature path

- **Every TRD physicochemistry column was the weighted mean over a subset of the repertoire.**
  `pchem_group` called `physchem_profile(group_by="locus")`, and the locus is derived **per row**
  from the gene name — so a TRD frame, which legitimately carries **TRAV** V genes because TRA and
  TRD share their V segments, split into two groups, and the loop over the tidy result let the last
  group win. Measured on a 10,000-clonotype synthetic TRD repertoire, 1,441 of whose clonotypes
  carried TRAV calls: all 30 TRD `pchem` columns moved on the fix, e.g. `all_charge`
  −0.00414 → +0.00124 and `center_volume` 97.792 → 97.857. Nothing raised, and the value was a
  plausible physchem mean throughout. The other six loci are unaffected (they carry one V naming
  scheme), and their agreement with the previous implementation is exact to 1.1e-12 absolute.
- **A docstring that cost a design decision.** `draw_pool` claimed the generator runs at 20,372
  sequences/s on IGH. Measured 2026-09-27 on one core of a 16-core M-series laptop,
  `generate(load_bundled("IGH"), 50_000, productive_only=True)` runs at **4,200 seq/s**, and
  identically at 1 and 16 polars threads — it is a per-sequence Python loop and does not thread.
  The 4.8x error is what made pool generation look like a 16 s stage when it is a 77 s one.

### Performance

Measured on a 16-core M-series laptop; the cohort shape is stated because a signature timing without
one has been wrong here before by 9x in the wrong direction.

- **One sample, seven loci, 10,000 clonotypes each: 0.732 s → 0.273 s (2.68x).** The composition
  groups were Python loops over every residue of every clonotype — 1.2M iterations and 208,863
  scalar `np.clip` calls per sample. `aa`, `kmer`, `spec`, `vus`/`jus`, `iso` and the QC fallback
  fractions are now `bincount`/`reduceat` over a 256-entry residue lookup table, and `pchem` replaces
  an explode-to-one-row-per-residue polars pipeline (150,000 rows per locus) with a `reduceat`
  against a 20x15 property matrix. Every value is unchanged except the TRD fix above.
- **A whole build, 64 samples of 5,000 receptors across seven loci: 300.8 s → 50.2 s (6.0x).** The
  build now runs **across samples** in worker processes: each worker draws its own repertoires from
  memory-mapped pools and returns one row of numbers, so nothing larger than a feature row crosses a
  process boundary — shipping drawn repertoires instead would be ~38 GB of pickling at the shipped
  size. Three faults were measured and fixed in that stage, in this order: pool generation fanned
  out over the **locus count** (7, on any machine) rather than the core count, 296.6 s; each of
  those workers then started a kernel with every core, 112 threads on 16 cores, 122.2 s after the
  first fix; and pool tasks were handed out in contiguous per-worker slices, which leaves workers
  idle behind IGH's pool, 51.1 s after the second. Featurisation of that shape is 2.5 s of the 50.2.
- The artifact is **bit-identical at every `n_jobs`**, because a sample is a pure function of its
  index: each one gets its own generator seeded from `(seed, locus, j)`. `POOL_CHUNKS` is a recorded
  constant, not the core count, for the same reason — a chunk count read off the machine would make
  an artifact built on 64 cores differ from the same build on 16. Pinned by
  `tests/python/test_signature_vectorized.py`.
- `vdjtools corpus --jobs/-j` (and `mir corpus -j`): worker **processes**, `0` = the whole
  allocation. The help text says which layer it reaches, because wiring a `--threads` flag to
  `n_jobs` is a mistake this repo has already shipped once.

### Added

- `vdjtools.signature.corpus` — build a synthetic corpus, winsorize it, fit the rotation and the
  scaling, read and write the artifact. `synthesize`, `fit`, `fit_locus`, `apply`, `Corpus`.
- `vdjtools.signature.features` — raw features and channels for one sample, a pure function of that
  sample plus the germline vocabulary.
- `vdjtools corpus` — the builder as a command, with `--smoke` for a minutes-long reduced build.
- New raw feature groups: **V usage**, **J usage**, **spectratype** (junction length per V gene) and
  **2-mer composition**, all clr/arcsine and folded into the same per-locus rotation. Raw widths run
  685 (TRD) to 3,744 (IGH), 13,483 features in total.
- New channels: `vsig:cov:<locus>:cstar` — the coverage the sample actually attained, **always**
  emitted, including when the diversity features are holes. It was previously computed, used to
  decide whether 28 columns were holes, and discarded: the only vsig quantity that was measured and
  thrown away. And `vsig:qc:-:winsor_frac`, what fraction of the row the corpus's bounds clamped.
- `support` on every raw feature declaration, as a closed vocabulary
  (`nonneg`/`nonpos`/`real`/`unit`), and winsorization **by percentile** with the side read from it.

### Removed

- **`vsig:pgen:*`** and the `pgen_q05` corpus constant. Per-clonotype Pgen was 96% of this command's
  runtime — 1,083 s of a 1,130 s cohort, essentially all of it IGH. Pgen remains a first-class
  `vdjtools.model` API; it is no longer a signature feature.
- `signature/blocks.py`, `assemble.py`, `presets.py`, `kmer.py`, `features/kmer_space.py` and the
  `vdjtools presets` command. The frozen TF-IDF + SVD k-mer space is subsumed by the corpus rotation.
- The reference-rescaling half of `signature/transform.py` (`reference_z`, `robust_loc_scale`,
  `magnitude_scale`, `DEFAULT_CLIP`). The transforms themselves are unchanged.
- **Tiers.** `core`/`standard`/`full` traded width for cost, and with a joint per-locus rotation
  every group must be computed before any component exists, so a tier could no longer skip work.
  `--components` replaces it and does the job better.
- CLI flags `--preset`, `--tier`, `--pgen-n-max`, `--cstar`, `--channels`, `--threads`. Each is
  pinned by a test that asserts a non-zero exit, because a flag that is accepted and silently
  ignored is the worst of the three outcomes.

### Changed

- **A corpus is required.** No default, because a silently chosen rotation makes two matrices look
  comparable when they are not.
- **The coverage target is a runtime argument** (`--cstar-target`, default: this cohort's own
  per-locus minimum attained coverage), not a constant in an artifact. A per-sample quantity in a
  corpus artifact is what produced the `tissue` defect above.
- `--components` takes a count or a variance fraction. The artifact stores the rotation up to the
  fitted count plus the **full eigenvalue spectrum**, so truncating downward later is exact and
  needs no refit, while asking for more raises and quotes what the spectrum reaches.
- `--jobs` is processes and says so. The old `--threads` reached `n_jobs`.
- `columns=` selects output rather than skipping work, for the joint-rotation reason above.
  Declining a whole locus still skips it.

### Measured while building this

- **A synthetic corpus must draw its depth.** At one fixed size, `depth:reads`, `depth:richness` and
  all five `pair:` log-ratios are identical in every sample, so their corpus spread is 0, they
  contribute nothing to the rotation, and the cross-locus block **cannot be fitted at all** — which
  fails quietly, as an artifact that simply omits it. Depths are now drawn log-uniformly over the
  measured real p05–p95 spread per locus (4.1x on TRB, 11.0x on IGH).
- **`numpy.random.Generator.zipf` is the wrong Zipf.** It samples integers *from* a Zipf
  distribution, which at `a = 1.5` has infinite mean, so normalising a draw gives one clone almost
  all the mass: repertoires collapsed to as few as 1 surviving clonotype and every
  coverage-standardised diversity feature became a hole. The spec is the Zipf law over **ranks**,
  `f_i ∝ i^-a`, which is a well-behaved rank-abundance curve: 336 of 500 clones surviving with 129
  singletons at 20 reads per clone.
- **Generation throughput**, 16-core M-series, one core: IGH 20,372 seq/s, TRB 23,113, TRA 34,506 —
  so a 10^7 pool is 8–14 min per locus, ~10 min for all seven at one process per locus.
- **Corpus build is byte-identical across thread counts**, verified `cmp` on the npz with
  `OMP_NUM_THREADS`/`POLARS_MAX_THREADS` at 1 against the default.

### Fixed: a functional V gene whose germline was one codon off scored Pgen exactly 0

`TRBV4-3*02`'s CDR3-region germline in both bundled human TRB models began `CTCTGCGCCAGC…` — one
whole framework codon (`CTC`, Leu) **upstream of Cys104**. OLGA's `human_T_beta` records its anchor as
267 against its own **287**-nt germline while `*01` is **284** nt at the same 267, so the defect is in
OLGA's tables and we reproduced it faithfully.

Reproducing it was not harmless, because `collapse_alleles` ranked a gene's representative germline by
**length** first: `*02`'s broken cut is 24 nt against `*01`'s correct 21, so the collapsed `TRBV4-3`
inherited the mis-anchored germline relabelled `*01`. Measured on 25,000 real human TRB clonotypes
(`isalgo/airr_control` `human.trb.ntvj`), `pgen_aa` was **exactly 0 for 698 of 698** TRBV4-3
junctions — **2.8%** of the set — with no error raised, against **1/595** zeros for the `TRBV4-1`
control. OLGA itself returns **1.063e-08** for `CASSQDLNTEAFF | TRBV4-3 | TRBJ1-1`: it marginalises
over the gene's alleles and `*01` carries the mass, so OLGA never returns 0 there. Fixing this
**restores** the exact-OLGA invariant rather than breaking it.

Two independent defences, because either alone leaves the other's failure reachable:

- **`model.io.repair_anchors`** makes the conserved anchor an invariant, on build **and on load**, so
  every already-shipped artifact is corrected on read with no regeneration. A V anchor steps by whole
  framework codons until the region starts at Cys104; a J is repaired only by undoing the known
  wrong-side slice (`full[anchor:]` where a J needs `full[:anchor+3]`), never by searching — a search
  finds an earlier in-frame Phe by coincidence and *shortens* the germline, which strands the allele's
  own deletion mass. Measured: on mouse `TRAJ19*01` a two-codon search cut arda's 30-nt germline to
  24. Repaired across the bundled set: `TRBV4-3*02` in `olga/TRB` and `learned/TRB` (anchor 267 → 270),
  and **9 of the 11** `learned` human TRA J alleles carrying framework downstream of Phe118 — 7 of
  those now match arda's germline **exactly**, a defect previously recorded as needing an all-loci
  regeneration.
- **`collapse_alleles` ranks the frame gate above length**, so a longer out-of-frame germline can
  never represent a gene. The gate is inert where no candidate passes it, which is load-bearing: an
  empty germline is out of frame too, and promoting it over a non-empty one is the `TRBV23/OR9-2`
  trap from the other side.

Consequence for the amino-acid scenario DP on the same 25,000 clonotypes: junctions the model could
not explain at all fell from **731 to 33**, 698 of them this one gene.

Not fixed, and not guessable from our reference: `TRAJ35*01` (human), the survivor of that TRA family.
arda's recorded anchor (25) points at a **Cys** codon, so neither the stored slice nor arda says where
Phe118 is; settling it needs IMGT arbitration rather than a search. Sizing is unchanged from 3.9.1:
**0** of VDJdb's 30,937 human TRA records use any of the family.

### Added: the amino-acid scenario API is usable from outside (#179)

`best_aa_scenarios` enumerates the plausible recombinations of an amino-acid junction with `len_v`,
`len_j` and D geometry. Three things stopped it being reachable as that primitive, all addressed:

- **Exported at `vdjtools.model`**, beside `infer_nt` and `best_scenario`, with
  `best_aa_scenarios_batch` and `gene_to_allele`.
- **Gene-level calls resolve** to a representative allele (`resolve_genes=True`, as
  `sc.paired_pgen`), because every real V/J call is gene-level and every one of them used to raise.
  A call naming no gene the model carries still raises and **names it**. The resolver now has one
  definition — `sc.pgen` carried a private copy, and #179's point was that every caller writes its
  own.
- **`best_aa_scenarios_batch`** returns a `polars.DataFrame`, one row per scenario, parallelized
  across sequences in native code. Measured on 24,907 real human TRB clonotypes at `k=8`:
  **8.02 s → 0.64 s** (12.5x), identical to the per-row loop and **bit-identical at any thread
  count**. A query the DP cannot explain contributes no rows, so an absent `row` is how a decline
  shows up — never a silently V/J-marginalised scenario.
- **Documented for what it is for**: the alternatives with their probabilities, D geometry, and naming
  a missing V/J by marginalising. It is *not* the better way to place a known boundary — measured on
  the same rows, top-1 places `v_end` exactly 59.8% and `j_start` 90.3% against germline alignment's
  73.3% / 97.7% at a twelfth of the cost.

### Changed: a corpus artifact is gated on the germline it was drawn from, not the library version

`Corpus.verify` compared `vdjtools.__version__`, which is both too strict — a patch release
invalidated every corpus — and wrong: it is read from installed distribution metadata, so a build
driven by `PYTHONPATH` against a different installed version records *that* version. The first
cluster artifacts recorded `3.6.0` from a 4.0.0 tree. The manifest now carries `models`, a per-locus
hash of the **cut segments the generator will actually draw from**, which a declared version string
cannot track: repairing `TRBV4-3*02` changed what every TRB pool contains while leaving
`olga:human_T_beta@2.0.0` identical.

### Measured: what the shipped corpora actually are

Two tables the docs now carry, both read off the artifacts rather than asserted
(:ref:`sig-component-table`, :ref:`sig-depth-sweep`).

**The two `--components` forms are not interchangeable.** Reaching 0.90 cumulative variance needs
**220** (TRG, `memory`) to **840** (IGH, `memory`) components on the statistics half, so a fraction
raises on the shipped `vsig` corpora rather than returning a narrower matrix — the rotation genuinely
stops at 128, and the message names the fraction reached. On the geometry half 0.90 needs **3 to 28**,
so a fraction is the natural knob there. IGH needs *more* components under `memory` than `naive` (840
against 676) while every other locus needs fewer.

**A corpus is depth-portable on the geometry half and not on the statistics half.** Twelve extra
corpora at N=1,000, drawn at the per-locus p05 / median / p95 of real `n_eff`, both halves, both
regimes. On `vsig` the centre moves **1.67 to 6.48 p05-robust-SD** from the shallow to the deep point
and per-PC scale falls to **0.285–0.904** of its shallow value; on `rsig` the centre moves
**0.03–0.08 SD** and per-PC scale stays within **±5%**. An `rsig` coordinate is a weighted mean of
fixed vectors, so depth changes its variance and not its value. The exception is `rsig` `memory`
bounds, which widen **4.6–7.8x** where `naive` stays at 0.94–0.98: bounds are percentiles, so clonal
expansion fattens the tails while leaving the centre and scale put.

**The builder is byte-reproducible.** Two smoke corpora built with different thread counts
(`POLARS_MAX_THREADS` 1 against the default) and different worker counts (`--jobs 3` against the
allocation) are `cmp`-identical. This is an acceptance criterion, not an aspiration: a corpus whose
values depend on the builder's machine cannot be compared with one built anywhere else.


## 3.18.1 — 2026-09-26

Guards for the class of bug 3.18.0 fixed, and one file that should never have been public.

### Added — `tests/python/test_knob_audit.py`

Every knob that moves a number must move only the columns it is documented to move. Turn the
knob, take the set of columns that moved, check it against the set the knob may touch. Seven loci
on purpose: a TRB-only fixture cannot see a reference that covers TRA and TRB and leaves five loci
of diversity as holes, which is the case that used to raise `KeyError`.

Mutation-checked rather than assumed — each of these fails the suite:

| mutation | caught by |
|---|---|
| a partial `cstar` dict borrows another locus's level | the TRA/TRB-reference test |
| `cstar_fallback_frac` hard-coded to 0 (the column exists but lies) | the blast-radius test |
| the coverage level leaks into `vsig:len` | two tests |

At cohort scale on 40 seven-locus samples: switching from the flat `C* = 0.20` to a measured
`C* = 0.1072` moves **992 cells, 18.14% of the matrix, across 29 columns**, and switching to a
TRA/TRB-only reference turns **800 of 1,120** `vsig:div` cells into declared holes with 200
`estimable` flags going to 0. Before 3.18.0 the first was silent and the second was a `KeyError`.

### Removed — `SOURCES.md` is no longer tracked

It is the one markdown file whose *purpose* is internal detail: cluster project paths, private
HuggingFace dataset names, local git-LFS checkout locations, and cohort composition down to batch
names and HLA schemas. It was already excluded from the sdist; it is now gitignored too, and the
references to it in `bundled.py`, `convert.py` and `SKILL.md` point at the per-artifact
`manifest.json` instead, which is where an outside reader can actually look.

Two CHANGELOG measurements were also reworded to name what was measured rather than where: a
cluster node instead of its hostname, "a bulk-RNA-seq AIRR cohort" instead of a project's store.

### Changed

`vsig`'s fallback bookkeeping is a counter, not a dict whose keys always equalled `reads`'.
`vdjtools signature --cstar none` now has a CLI test.

## 3.18.0 — 2026-09-26

The coverage level the diversity block rests on is now **visible in the vector**. It was not, and
that is the whole of `ISSUES.md` item 6 (3): `DEFAULT_CSTAR = 0.20` is one number stretched over
seven loci that attain very different coverage, and a `vsig:div` block computed at it is fully
populated and entirely plausible. Nothing said it was a guess.

### Added — `vsig:qc:-:cstar_fallback_frac`

The share of a sample's **present** loci whose coverage level came from the flat fallback rather
than from a measured constant. One column, in every tier, in the `qc` channel beside the V/J
fallback fractions — for the same reason that block exists: it is the number that tells a
collaborator whether your vector is comparable to theirs.

The size of what it declares, 40 synthetic TRB samples of 400–4,000 clonotypes, median absolute
shift against the bulk-RNA-seq measured level (`C* = 0.1072`), in `log10` Hill units:

| coverage level | `1D_c` | `0D_c` | `2D_c` | `clonality` |
|---|---:|---:|---:|---:|
| flat fallback `C* = 0.20` | 0.2684 | 0.2765 | 0.2575 | 0.4575 |
| amplicon `C* = 0.4080` | 0.5755 | 0.6030 | 0.5432 | 1.0329 |

0.2684 in `log10` is a **factor of 1.86** on the Shannon diversity, and against the amplicon level
a factor of 3.8. Every sample in both columns is estimable (40/40), so none of this shows up as a
hole; two matrices differing by that much join cleanly and look the same.

Tier widths move by one: `core` 152 → **153**, `standard` 688 → **689**, `full` 1403 → **1404**;
the `vsig` half 160 → **161**. Only the `nuisance` preset gains it (73 → 74) — it describes the
run, not the donor, so `classify` (615), `transfer` (550), `compact` (86), `bcell` (271),
`geometry` (514) and `statistics` (101) are unchanged. The shipped scale references cover 1,403
columns and do not scale the new one; it passes through, which is right for a provenance fraction.

### Changed — a locus with no established level is a hole, not a borrowed number

`cstar` takes `None` and takes a **partial** dict. A locus the dict does not cover establishes no
level: its `vsig:div:*` columns are `nan` and `vsig:mask:<locus>:estimable` is `0`, with one
warning naming the loci. Standardising IGH's Hill numbers to a level measured on TRB is not a
measurement, and a plausible number in the wrong units cannot be detected by the caller while a
`nan` can.

The refusal and the fallback are **disjoint**, so the two signals add up: a locus with no level is
counted by `estimable`, a locus with a borrowed one by `cstar_fallback_frac`. `cstar=None` reports
`0.0`, not `1.0` — it did not fall back, it declined.

`vdjtools signature --cstar none` reaches the same path from the command line, and the stderr note
says which of the two the run used.

### Fixed — a partial `cstar` dict raised `KeyError` from inside the loop

`level = cstar[locus]` indexed the dict directly, so a reference covering TRA and TRB only — which
is exactly what the shipped amplicon reference is — died on the first B-cell locus instead of
reporting a hole. `mir.signature.signature()` worked around it by completing the dict with a
deliberately unreachable `C* = 1.0` so the estimability check would fail; that workaround is gone
in mirpy 3.20.0 because the library now owns the semantics, and the wasted `estimate_d` call per
uncovered locus goes with it.

### Docs

`docs/signature.rst` gains the provenance table above, in the section that already explains why
`C*` is load-bearing. The preset table's `compact` (152 → 86) and `bcell` (286 → 271) widths were
stale before this release and are now the live values.

## 3.17.1 — 2026-09-26

The Pgen block you cannot decline. 3.16.0 let a caller skip `vsig:pgen` entirely, which is 96.6% of
a `standard`-tier seven-locus sample — but `--preset classify` and `--preset transfer`, the two
`recommended` presets, both keep all seven loci of it and paid the full bill. The native D DP now
costs 2.4-3.6x less on the D-bearing loci, and **every Pgen value it returns is bit-for-bit what
3.16.0 returned.**

### Changed — the D germline is threaded once per 5' trim, not once per 3' trim

`pgen_aa_vdj` enumerated `(D allele, ndel5, ndel3)` states and re-threaded the surviving germline
through the 25-state codon DP from scratch for each one. The emission at `idx3` is a **prefix** of
the emission at `idx3 - 1` — the trims of one `(D, 5' cut)` are nested — so one forward pass now
covers every 3' trim of that chain and the `len(D)` factor leaves the cost. The argmax mirror
`best_vdj` has always done this; the sum never did.

Inner nucleotide steps at a 20-aa junction, counted off the bundled models:

| locus | D states, p > 0 | (D, 5' trim) chains | steps before | steps after | |
|---|---:|---:|---:|---:|---|
| IGH | 9,212 | 767 | 4,736,568 | 588,010 | **8.1x** |
| TRD | 455 | 45 | 139,568 | 22,095 | 6.3x |
| TRB | 297 | 38 | 107,087 | 21,147 | 5.1x |

Wall clock, 250 junctions per locus drawn from the bundled `olga` model, 16-core M-series:

| locus | mean junction | before, all cores | after | | before, 1 thread | after | |
|---|---:|---:|---:|---|---:|---:|---|
| IGH | 20.9 aa | 2,557 µs | 703 µs | **3.6x** | 30,106 µs | 7,859 µs | **3.8x** |
| TRD | 18.8 aa | 155.3 µs | 64.8 µs | 2.4x | 1,400 µs | 557 µs | 2.5x |
| TRB | 15.1 aa | 79.6 µs | 27.6 µs | 2.9x | 666 µs | 227 µs | 2.9x |
| TRA | 14.2 aa | 113.0 µs | 108.5 µs | — | 923 µs | 912 µs | — |
| TRG | 13.5 aa | 13.0 µs | 13.0 µs | — | 109 µs | 110 µs | — |
| IGK | 11.4 aa | 9.1 µs | 8.8 µs | — | 67.0 µs | 68.6 µs | — |
| IGL | 12.9 aa | 8.4 µs | 8.1 µs | — | 55.4 µs | 57.2 µs | — |

The four VJ loci have no D and are the control: they do not move, which is what says the change is
where it claims to be.

End to end, which is the number the presets care about — **120 synthetic seven-locus samples** at
the real cohort shape (median 1,388 clonotypes, 39x p05-p95 spread, 30% IGH / 27% TRB / 17% IGK /
13% TRA / 10% IGL / 1% TRG / 0.3% TRD, 260,158 clonotypes in total), `--preset classify`, the
`vsig` half, `threads=0`:

| | before | after | |
|---|---:|---:|---|
| 120 samples, total | 198.5 s | 64.6 s | **3.1x** |
| per sample | 1,654 ms | 539 ms | |

A TRB-only benchmark cannot see any of this: IGH is ~96% of the bill and the other six loci
together are the rest.

Conditioning Pgen on the observed V/J is still not a lever, re-measured on this release: IGH 1.1x,
TRB 1.1x, TRD 0.9x. It remains a large win only on the four D-less loci, which together cost under
2% — and it would shift `log10 Pgen` by about -1 decade on IGH against a frozen reference. Slower
where it matters and wrong everywhere.

`combine_tm` is untouched and is now the bottleneck — ~437k calls per 20-aa IGH junction against
588k threading steps — which is why the wall clock gains 3.6x where the steps gain 8.1x. Hoisting
its per-boundary reduction is worth roughly another 2x and is **not** taken here: it regroups the
sum, and that forfeits the guarantee below. It is noted at the call site with its cost.

### The output does not move, and that is checked rather than asserted

The rewrite is arranged in two passes: the walk records each join, then a second pass adds them into
the total **in the original `(ndel3, position)` order**. Summing straight out of the walk would be
correct to ~1e-16 and would still be wrong to ship, because `vsig:pgen:*:frac_atypical` is measured
against a `pgen_q05` reference frozen outside the block — a last-bit move there is a silent data
bug that lands in a plausible range, not a rounding detail.

- **1,708 of 1,708 Pgen values bit-for-bit identical** on one platform, 3.16.0 against 3.17.1,
  across all 7 loci x both `from_olga` (uncollapsed) and `load_bundled` (collapsed, the default) x
  `pgen_nt`, `pgen_aa`, `pgen_aa(mismatches=1)`, V/J-agnostic and V/J-restricted. Collapsing
  changes the D-state cardinality the DP walks, so both settings are covered. Same platform is the
  scope claimed and the scope that matters: upgrading does not move your numbers.
- **7-locus OLGA concordance reproduces `r(log10 Pgen) = 1.00000`** on nt and aa. IGH agrees with
  OLGA to a max relative error of **1.0e-13 (nt) / 7.9e-14 (aa)**.

Across *compilers* it is a different question with a pre-existing answer: float64 Pgen has never
been bit-portable. The same bundled model and the same junction return a **1-ULP** (~2e-16)
different value under Linux/GCC than under macOS/clang, because the sum contracts differently — no
code change involved, and true of every release before this one.

### Added — `tests/python/fixtures/pgen_golden.json`, the first stored Pgen reference

276 values over IGH/TRB/TRD as hex float64, compared at a max relative error of **1e-9** — the bar
`appendix/compare_models.py` already calls "EXACT". It exists because the oracle tests measurably
do not catch this class of bug: stopping the shared walk one nucleotide short — each chain loses its
longest D — moves IGH Pgen by a **median relative error of 7.9e-7** over 12 bundled-model junctions,
*under* the `rtol=1e-6` that OLGA comparisons here use, even though the worst junction moves 4.3%.
All three new IGH oracle tests passed against that mutation; the frozen reference caught it on every
value — 800x tighter than the signal, and seven orders looser than cross-compiler drift. Regenerate
it only when a germline or a bundled model changes.

`test_pgen_d_prefix.py` also closes a gap this change made untenable: **IGH had no OLGA comparison
anywhere in the suite**, though it is the 9,212-state locus and 96% of the signature bill. It now
has one at `rtol=1e-9`, nt and aa, on a length-spread draw rather than the shortest junctions.

### Unchanged on purpose

Three other sites carry the same nested-prefix structure and each now says at its own definition why
it was left: `pgen_aa_vdj_dd` (no bundled model has `P(n_D=2) > 0`, and its D1 loop accumulates into
a shared message array so order preservation is a different problem), `d_middle` (its match loop
exits on the first mismatching nucleotide, so the redundancy is not its cost — the V x J enumeration
around it is), and `accum_vdj` (soft counts, pinned at `atol=1e-12` rather than bitwise, and it runs
only during model regeneration).

### Note on 3.17.0

`v3.17.0` was tagged with this same change and **never published** — its pre-publish CI job went red
because the new frozen reference was compared with `==`, asserting bit-identity across compilers
rather than across versions. Nothing reached PyPI under that version; 3.17.1 is the first release of
this work.

## 3.17.0 — 2026-09-26

### Changed — the D-germline walk is shared across 3' trims

`pgen_aa_vdj` enumerated `(D allele, ndel5, ndel3)` states and re-threaded the surviving germline
through the 25-state codon DP from scratch for each one. The emission at `idx3` is a **prefix** of
the emission at `idx3-1`, so the trims of one `(D, 5' cut)` are nested: one forward pass now covers
every 3' trim of that chain and the `len(D)` factor leaves the cost. `best_vdj`, the argmax mirror,
had always done this; the sum never did.

| locus | measured |
|---|---|
| IGH | 9,212 live D states are 767 prefix chains; inner nucleotide steps at a 20-aa junction 4,736,568 -> 588,010 (8.1x); 250 junctions 2,557 -> 703 microseconds each (**3.6x**) |
| TRB | 2.9x |
| TRD | 2.4x |
| the four D-less loci | unmoved |

End to end on 120 synthetic seven-locus samples, `--preset classify`: 198.5 -> 64.6 s, 1,654 -> 539
ms/sample (3.1x).

**Output is bitwise unchanged, by construction and by measurement.** The walk records each join and
a second pass adds them in the original `(ndel3, position)` order; summing straight out of the walk
would reassociate, and a last-bit move is a silent data bug where a Pgen is compared against a
quantile frozen outside the block. Verified on **1,708 of 1,708** float64 values across 7 loci x
uncollapsed and collapsed models x `pgen_nt` / `pgen_aa` / `pgen_aa(mismatches=1)`, agnostic and
V/J-restricted. The seven-locus OLGA concordance still reproduces `r(log10 Pgen) = 1.00000`, IGH at
max relative error 1.0e-13.

Adds `tests/python/fixtures/pgen_golden.json`, the repo's first stored Pgen reference, compared with
`==` rather than a tolerance — because a tolerance measurably does not catch this bug class:
stopping the walk one nucleotide short moves IGH Pgen by a median **7.9e-7** relative, under the
`rtol=1e-6` every OLGA test here uses, while the worst junction moves **4.3%**. The same fixture
closes a longer-standing gap: IGH, the 9,212-state locus, had no OLGA comparison anywhere in the
suite.

Three sites carrying the same structure are deliberately left alone and each now says why at its own
definition: `pgen_aa_vdj_dd` (no bundled model has `P(n_D=2) > 0`), `d_middle` (its match loop exits
on the first mismatch, so the V x J enumeration around it is the cost) and `accum_vdj` (soft counts,
pinned at `atol=1e-12`).

## 3.16.0 — 2026-09-26

The signature hot path, measured rather than assumed. Headline: **`vsig:pgen` is 96.6% of a
`standard`-tier seven-locus sample**, essentially all of it IGH — and you can now decline it.

### Added — `columns=` skips the work instead of dropping the result

`vsig(columns=...)` and `vsig_cohort(columns=...)` take an explicit subset, intersected with the
tier and returned in layout order. **A block whose columns are all excluded does not run.** That
distinction is the whole feature: a select-and-drop would look identical in the output and save
nothing.

| what | before | after | |
|---|---:|---:|---|
| `standard` tier, one sample, no `vsig:pgen:` | 6,234 ms | 211 ms | **29.6x** |
| `--preset nuisance` (`full` tier, zero pgen columns) | 1,283 ms | 41 ms | **31.0x** |

Every column kept is bit-identical to the same column from a full run — pinned by
`test_signature_columns.py`, which booby-traps each gated block in turn rather than timing it (a
stopwatch cannot tell a skipped block from a fast one).

`vdjtools signature --preset` now passes the preset's columns into the computation. It was a
display filter; `nuisance` computed all seven loci of Pgen and threw every one of them away.

### Added — `pgen_n_max` / `--pgen-n-max`

How many junctions per locus the Pgen block measures. Cost is linear in it. **Not free**: mean and
sd are unbiased at any *n*, but `pgen:*:frac_atypical` is compared against a frozen `pgen_q05`
reference drawn at 2,000, and a different draw is a different quantity. The default is unchanged.

### Fixed — the documented low-memory path did not work

`vsig_cohort` and mirpy's `signature_cohort` both advertised a zero-argument callable per sample as
the `O(n_jobs)`-memory path, but only mirpy's `_one_rsig` ever resolved one. Passing a callable
died on `AttributeError: 'function' object has no attribute 'columns'`, so the documented path was
unusable rather than merely slow. `signature.cohort.resolve_sample` is now the single definition,
used by all three entry points.

### Fixed — `signature --describe` and `--channels` described the other tool's output too

Since the 3.18.0 mirpy split this command emits `vsig` only, but `--describe` still printed all 688
columns and `--channels` all 20 channels. The one output whose entire job is "the exact columns you
will get" was naming 528 columns you will not get. Both are filtered to this half; `mir signature`
is the mirror image.

### Changed — the Pgen cost is documented, and three proposed fixes are recorded as dead

`pgen_block`'s docstring carried a diagnosis that measurement does not support: that V and J are
marginalised for an API reason, read as though conditioning on the observed call were the blocked
speedup. It is not. 250 junctions per locus, bundled models, `threads=0`:

| locus | marginalised µs/junction | conditioned | speedup | D states, p > 0 |
|---|---:|---:|---:|---:|
| IGH | 2,675 | 2,767 | **1.0x** | **9,212** |
| TRD | 146 | 135 | 1.1x | 455 |
| TRB | 85 | 74 | 1.1x | 297 |
| TRA | 108 | 4 | 26.9x | 0 |
| TRG | 16 | 4 | 4.0x | 0 |
| IGK | 10 | 2 | 3.9x | 0 |
| IGL | 10 | 2 | 4.3x | 0 |

Conditioning is a large win only on the four D-less loci, which together are under 2% of the bill.
IGH costs 31x TRB because it has 31x the `(D allele, ndel5, ndel3)` states — 9,212 against 297,
matching the cost ratio to the digit. It would also shift `log10 Pgen` by **-1.07 decades on IGH**
against a frozen reference, so it is slower *and* wrong.

Also measured and recorded as dead ends, so they are not re-derived: collapsing germline-identical
D states (9,212 states -> 8,910 distinct emissions, **1.0x**), and memoising per unique junction
(IGH deduplicates **1.0x** on real cohorts). A future win has to come from the native D DP, holding
the exact-Pgen invariant.

### Changed — the block that dominates once Pgen is declined, 4-17x faster

With the Pgen block gone, `div_block` is **84% of everything that remains**, and essentially all
of it was `stats.inext._rtd_moment` — a Python loop over the abundance *spectrum*, i.e. once per
**distinct clone size** (~200 on a real locus), called **28 times per sample**: orders q=0 and
q=1, floor and ceil of a non-integer `m`, seven loci. The loop body was a small numpy op, so the
cost was round trips, not arithmetic.

Now one flat pass: the ragged `k` ranges are flattened with `repeat`/`cumsum` and regrouped with
`add.reduceat`, so each distinct count's terms are still summed among themselves before the
groups are combined.

| `div_block`, seven loci | before | after | |
|---|---:|---:|---|
| 415-clonotype sample | 58.3 ms | 3.4 ms | **17.1x** |
| 4,863-clonotype sample | 286.9 ms | 70.8 ms | **4.1x** |

The flat buffer is bounded (`_FLAT_TERMS`) and processed in chunks: unchunked, a 1.5M-read
single-locus amplicon sample regressed ~18% against the old loop, because there the groups were
already large enough that the loop was doing big numpy ops. Chunked it is a wash on that shape
and 4-17x on ordinary repertoire samples.

Worst relative difference against the loop it replaces: **1.5e-15** over 24
(depth, order, size) combinations, and **1.0e-14** across every `vsig:div:*` column on real
samples — floating-point reassociation, nothing more. `test_inext_vectorised.py` pins the
equality rather than the time, because a speedup that moves a diversity value is worthless when
the value is compared against a frozen scale reference.

### Fixed — a pool worker asked for the whole machine's Pgen threads as well

`vsig_cohort(n_jobs=N)` passes no `threads`, so each of the N workers asked for
`hardware_concurrency - 2`: **224 Pgen threads on a 16-core box at `n_jobs=16`**. The pool and
the kernel each sized themselves off the machine with nothing reconciling them. `pgen_block` now
takes one thread when it is inside a pool worker (`signature.cohort.in_pool_worker`).

### Added — `signature` says that its output is unstandardised

`mir signature` standardises its half against a frozen scale reference. This command cannot: the
reference ships in mirpy and the dependency runs the other way. So a plain join of the two CLI
outputs is a **mixed-scale** matrix, and `vsig:pgen:*:frac_atypical` is `nan` without the
reference's `pgen_q05`. Measured on one synthetic 600-clonotype TRB sample at `tier="core"`:
`vsig:div:TRB:1D_c` reads 1.9502 standardised against 2.2495 raw.

That gap is better said than discovered, so the command now says it on stderr and points at
`mir.signature.signature_cohort()`, which applies the reference to both halves in one call.

### Changed — `_bundled_model`'s cost, corrected by a factor of 5-10

The docstring claimed 0.4-1.8 s per model, and that number was the stated justification for the
cache. Measured: **51-268 ms per locus, 1.12 s for all seven** (TRG 51, TRB 83, IGL 100, TRD 119,
TRA 239, IGK 259, IGH 268). The cache is still worth having and saves ~1.1 s per process — but it
is per *process*, so every spawned worker pays it again.

## 3.15.0 — 2026-09-25

### Added — one process pool, in one place

`vdjtools.signature.cohort` holds `slices()` and `parallel_rows()`: contiguous chunks, one task
per worker, `spawn` (polars cannot be combined with `fork`), and **no serial fallback** — a pool
that cannot start raises with the two fixes named. Three callers share it: `vsig_cohort` here,
`rsig_cohort` and `signature_cohort` in mirpy. Previously only mirpy had a pool and vdjtools had
none.

`vsig_cohort` gains `n_jobs` (default `1`, `0` = every core this process may use) and accepts a
zero-argument callable per sample, which defers the read into the worker and keeps peak memory at
`O(n_jobs)` samples rather than the whole cohort.

### Changed — the `signature` help no longer misdescribes mirpy

It said `mir signature --preset classify ...` "emits both halves as one vector". As of mirpy
3.18.0 it emits the geometry half only — one tool per half, joined on `sample_id`.

## 3.14.2 — 2026-09-25

### Fixed — pools were sized off the machine's cores, not this process's

`os.cpu_count()` reports the machine. It is not what a process is allowed to use, and everywhere
repertoire analysis actually runs, the two differ. Measured on a SLURM cluster node allocated
with `srun -c 8`:

| | |
|---|---:|
| `os.cpu_count()` | **40** |
| `len(os.sched_getaffinity(0))` | **8** |
| `os.process_cpu_count()` | 8 |

Every pool here sized itself off the first number, so an eight-core allocation started up to 40
workers. That is not merely wasteful: each worker pays a fresh interpreter and a fresh load of the
frozen artifacts, so on a 32 GB box it is several GB of overhead competing for cores that do not
exist, and it can take the machine out of memory on a cohort that would otherwise have fitted five
times over.

New `vdjtools.cores.available_cores()` takes the **smallest** of the affinity mask
(`os.process_cpu_count()` on 3.13+, else `os.sched_getaffinity`), the cgroup CFS quota, and
`os.cpu_count()`. Affinity and quota constrain independently and a box can carry both:
`docker run --cpuset-cpus` shows in the affinity mask, `--cpus=N` does **not** — it is a bandwidth
quota, invisible to every API except the cgroup file. Kubernetes CPU limits are the same
mechanism, which is why the cgroup read is there rather than trusting affinity alone.

Used at every pool-sizing site: `model/score.py` (`_pgen_nt_many`), `model/data.py` (`build_all`)
and `io/batch.py` (`map_samples`, which had been leaving `ThreadPoolExecutor` to its own
`min(32, os.cpu_count() + 4)` default — 36 threads for 8 cores on that node). A test walks the
package and fails on any new bare `os.cpu_count()`, so a future site cannot quietly re-create it.

## 3.14.1 — 2026-09-25

### Fixed — `on_duplicate` now reaches the signature path and both CLI commands

3.14.0's error message tells the caller to `pass on_duplicate="sum"`, and there was no way to do
it: `sanitise` took the argument, `vsig` did not, and neither `vdjtools diversity` nor
`vdjtools signature` exposed a flag. The only route through was to collapse every locus frame by
hand before calling, which is unworkable over a cohort. `vsig(..., on_duplicate=…)` now forwards
it (and so does `vsig_cohort`), and both commands take `--on-duplicate error|sum`. A message that
names a fix the caller cannot apply is barely better than the silence it replaced.

## 3.14.0 — 2026-09-25

### Fixed — a duplicated amino-acid clonotype key is now a hard failure

A frame with no `junction_nt` that repeats `(junction_aa, v_call, j_call, c_call)` cannot say
whether those rows are two nucleotide clonotypes encoding one peptide, or one clonotype the export
duplicated. Richness, clonality, Shannon and top-clone fraction all differ between the two
readings, and until now the library counted rows and picked the first reading in silence: the four
rows below returned `observed_diversity = 4`, with no warning, from `diversity_stats`,
`diversity_cohort` and the whole signature path.

```
CASSLGQGAYEQYF  TRBV5-1*01  TRBJ2-7*01  TRBC2  7
CASSLGQGAYEQYF  TRBV5-1*01  TRBJ2-7*01  TRBC2  3
CASSPRTGELFF    TRBV7-9*01  TRBJ2-2*01  TRBC2  5
CASSQDRGNTIYF   TRBV4-1*01  TRBJ1-3*01  TRBC2  2
```

**Measured cost of the silence.** Two exports of the same samples, one collapsed to the
amino-acid key and one not, went through the same estimators without complaint. Richness differed
by **1.5% overall and 5.0% in IGK**, and the affected samples then sat **3-5 robust-SD from the
rest of the cohort** on exactly the diversity and clonality columns — read as biology until the two
exports were diffed.

New in `vdjtools.io.schema`, and re-exported from `vdjtools.io`:

- `assert_resolvable(df, name=…)` — the gate. Raises, naming the sample, the offending keys, and
  both legitimate resolutions.
- `duplicate_keys(df)` — the repeated keys and how many rows each covers, worst first.
- `collapse_duplicates(df)` — the deliberate sum, recomputing `frequency`.
- `has_nt_resolution(df)` — whether `junction_nt` can tell two rows apart.
- `resolve_duplicates(df, on_duplicate)` — the policy dispatcher.

`diversity_stats`, `diversity_cohort` and `vdjtools.signature.blocks.sanitise` all gained
`on_duplicate="error"` (default) / `"sum"`. **The default is `error`**: the frame cannot answer the
question, so the library must not answer it either.

Two things deliberately do not trip the gate. A frame carrying `junction_nt` is never rejected —
the duplicates are then real and the key that resolves them is present. And an all-null
`junction_nt` column does **not** count as resolution: a schema-conformant frame always has the
column and an amino-acid-collapsed export fills it entirely with nulls, so presence alone would
have let the filed case straight through.

The cohort check stays lazy. It is one streamed `group_by` over
`(sample_id, junction_aa, v_call, j_call, c_call)`, filtered to the offending keys before
collection, and skipped entirely on a cohort that carries `junction_nt`.

### Changed — `_PREP_CACHE` is bounded

`model/viterbi.py`'s prepared-model cache is keyed on `id(model)` and already stored and verified
the model reference, so it never returned another model's entry. It was unbounded, which leaked one
pinned `Model` per model ever prepared. Now capped at 8. The bound is load-bearing twice: pinning
the model is also what keeps the `id()` key honest, since an id can only be reused once its owner
is collected — so eviction has to drop the id and the model together.

### Changed — the seqtree floor is 1.0.0

`seqtree>=0.6.1` predated seqtree's semver guarantee, which starts at 1.0.0, so the floor did not
actually bound what pip could resolve. Both the full suite (1,189 tests) and the docs build are
green against seqtree 1.0.0, vdjmatch 0.3.1 and arda-mapper 2.23.0; only the seqtree floor moves,
because the vdjmatch and arda floors each pin a specific fix and raising them would exclude working
versions for no stated reason.

### Audited — every pool and every cache, with measurements

Recorded in `CLAUDE.md` under "Pools and caches: the standing audit". All three thread pools wrap a
call that releases the GIL and all three are genuinely parallel:

| site | measured |
|---|---|
| `map_samples` | 2.24x / 3.49x / 3.85x at 2 / 4 / 8 workers (64 samples x 10k clonotypes) |
| `_pgen_nt_many` | 1.92x / 3.68x / 6.92x / 11.4x at 2 / 4 / 8 / 16 threads (256 TRB junctions) |
| `build_all` | structural — arda is a subprocess, `native.estep_batch` releases the GIL |

`map_samples` tops out near 4x rather than 16x because polars already multithreads each read; past
four workers the pool competes with polars' own thread pool for the same cores.

`tests/python/test_thread_scaling.py` pins the `_pgen_nt_many` claim so it fails loudly if the GIL
guard is ever dropped. It runs on **TRB** on purpose: a V/J-marginalized nt Pgen costs 0.01 ms on
TRG and 0.02 ms on TRA against 36.88 ms on TRB and 40.70 ms on TRD, so the same test written on a
VJ chain measures pool startup and reports a 0.83x *slowdown*.

Every `lru_cache` in the package loads a frozen artifact keyed on that artifact's identity, and
every one is bounded; none memoises a computation whose inputs are not fully in the key.

## 3.13.0 — 2026-09-24

### Added — the channel vocabulary

A signature column has always been named `<sig>:<block>:<locus>:<feature>`, and the second field
has always been the unit of interpretation. It had no name and no API, so a caller who wanted "the
diversity columns" reconstructed the grouping by string-splitting, and a finding came out as
"column 412 moved" rather than as a sentence.

The field is now called the **channel** and is addressable:

- `CHANNELS` — the vocabulary, 20 entries, each channel to the one quantity it measures.
- `channel(column)` — the channel a column belongs to.
- `channels(tier, sig, columns=…, per_locus=…)` — channel to column indices, over any column list,
  optionally keyed per locus. **Disjoint and exhaustive**, so per-channel shares sum over the whole
  vector with nothing left over; a test asserts the partition at every tier.
- `channel_table(tier)` — one row per channel: width, loci, attributability, what it measures.
- `vdjtools signature --channels` — the same table from the command line, reading no input.

A channel key carries its half (`vsig:div`, not `div`): `depth` and `div` are declared on **both**
halves and are different measurements of the same idea, not duplicates. `attributable` — whether
"which clonotypes drive this" is a well-posed question — is read off the `Block` where it was
declared at build time, never inferred from a name.

The map is what `mir.signature.channel_spec` feeds to `mir.explain.channel_report`, which is how a
model's "it separates the groups" becomes "IGH diversity and isotype composition carry it".

### Changed — the README leads with what most people run

Repertoire analytics now opens the README, with a **Where to start** table routing to the docs
page for each task; the recombination model engine moved below the Python API, and the portable
signature became its own clearly-marked **extended** section linking to
[Signature](https://docs.isalgo.dev/vdjtools/signature.html) and the new
[Channels](https://docs.isalgo.dev/vdjtools/channels.html) page.

### Fixed — the docs site showed the wrong version

`docs/conf.py` hard-coded `3.6.1`, so every page rendered that in its header through six minor
releases. It now reads `vdjtools.__version__`.

### Fixed — the coverage constants quoted in the docs

The `cstar` table for TRA/TRB compared 0.545/0.408 against 0.129/0.126, quoted from a pre-3.11
fit. The references that actually ship carry 0.1230/0.1072 for bulk RNA-seq, so the amplicon-vs-
RNA-seq gap on TRB is **3.8x**, not 3.2x.

## 3.12.1 — 2026-09-24

### Fixed — a k-mer space no longer loads through `pickle`

`save_kmer_spaces` stored the locus list and the V-gene list as `dtype=object`, so
`load_kmer_spaces` had to pass `allow_pickle=True` — which makes **reading a space someone sent
you arbitrary code execution**. These files exist to travel; that is their entire purpose, so the
loader has to be safe on a file it did not write. Strings are now fixed-width unicode (`<U`) and
the loader passes `allow_pickle=False`.

A space written before this release carries object arrays and will now raise on load, with numpy's
own message naming the cause. Re-save it with `save_kmer_spaces` to convert.

### Changed — spaces are stored in float32, halving them

The IDF vector and the rotation are stored as float32 rather than float64. Measured across all
seven loci, projecting a repertoire through the float32 basis instead of the float64 one moves the
result by at most **4.9e-9 relative** — nine orders of magnitude below the noise in any repertoire
measurement — and it takes the fitted artifact from **14.2 MiB to 7.0 MiB**, which is the
difference between shipping it inside a wheel and not shipping it at all.

## 3.12.0 — 2026-09-24

### Changed — the docs tell you where the signature is

The signature pages were written, committed and live, and nobody could find them: the left sidebar
was a flat wall of page titles in source order, so `signature` sat in it indistinguishable from
everything else. `index.rst` now carries captioned toctrees — **Start here**, **Repertoire
signatures**, **Models and single cell**, **Worked examples**, **Reference** — and the theme is
configured to render that structure rather than flatten it, matching `seqtree` and `mhcmatch`.
Adds the `site-nav` sidebar template and the card CSS those two already use.

### Added — the end-to-end recipe, and a notebooks page

`docs/signature.rst` now answers the question people actually ask: a directory of per-sample AIRR
TSVs plus your own metadata sheet, to one table you can join. `sample_id` is the file name up to
the first dot, so naming files after the key your metadata already uses means the join needs no
mapping table.

`docs/notebooks.rst` links the example notebooks, which the docs did not reference at all.

### Changed — no worktrees

`CLAUDE.md` now says branches only. Work committed under a gitignored `.claude/worktrees/` is
invisible from the repo root and effectively unbacked, and a branch held by a worktree cannot be
checked out in the main repo.

## 3.11.0 — 2026-08-20

### Added — `vdjtools correct-vj`, and a documented pre-processing surface

`correct_vj_usage`/`apply_vj_correction` existed in the Python API and had no CLI face, so the
one operation you most want to run over a directory of samples was the one you had to write a
script for. `correct-vj` takes the samples and their batch labels, writes the corrected usage
table, and — with `--outdir` — rewrites each clonotype table:

```bash
vdjtools correct-vj s1.tsv s2.tsv s3.tsv s4.tsv -b A,A,B,B \
    --transform sigmoid --usage-out usage.tsv --outdir corrected/
```

`-b` takes a comma-separated list in file order, or a TSV with `sample_id`/`batch` columns
matched on the file stem. Two transforms: `location`, the location term of ComBat on usage
log-probabilities, and `sigmoid`, the σ-standardised grand-mean-preserving z-score of Vlasova
et al. 2026, which corrects a batch that is merely *noisier* in a gene — something a location
adjustment cannot do.

The command's help carries the warning the method needs: **name the technical variable.** Point
`--batches` at a primer mix, a run, an extraction protocol. Pointed at a study identifier that is
collinear with the biology being measured, it removes the effect along with the batch and nothing
reports that it happened.

### Added — `docs/preprocessing.rst`

One page for the whole pre-processing surface, which was previously discoverable only by reading
the API reference: format conversion, the three filtering axes from 3.10.0, how `frequency` is
handled at every step and when it is renormalised, error correction, downsampling, V/J-usage batch
correction, pooling and joining, the CLI equivalents, and a short section on how mirpy differs
(it applies the productive filter on every read and will not let you turn it off).

`README.md`, `docs/index.rst` and `docs/usage.rst` now point at it, and the three filtering axes
are spelled out in the README's CLI and Python quick-reference blocks rather than left implicit.

### Changed — examples updated to the 3.10.0 names

`cdr_features.py`, `overlap_similarity.py` and `preprocess.py` used `filter_functional(keep=
"coding")`. They now use `filter_productive`, and `preprocess.py` demonstrates all three axes plus
`correct_vj_usage`.

### Fixed — a private dataset now says so

`vdjtools model build` trains on `isalgo/airr_model_read`, which is private, so outside the lab it
stopped at a bare `huggingface_hub` 401 that named neither the repo nor what to do. It now raises a
`PermissionError` explaining that this is expected, that the bundled models are the versioned
reference and need no rebuild, and it re-raises anything that is not an auth failure so real errors
are not swallowed. The command's help says LAB-ONLY.

## 3.10.0 — 2026-08-20

### Added — three filtering axes, named after the standards that define them

One English word had been doing three jobs. They are now three predicates with three names:

| axis | question | standard | function |
|---|---|---|---|
| parseable | is `junction_aa` readable at all? | — | `assert_parseable` (raises) |
| productive | does the rearrangement encode a chain? | [AIRR](https://docs.airr-community.org/en/latest/datarep/rearrangements.html) | `filter_productive` |
| functional | is the germline gene real? | [IMGT](https://www.imgt.org/IMGTindex/functionality.php) F/ORF/P | `filter_functional_genes` |

AIRR's `productive` is a property of the **rearrangement**; IMGT's *functional* is a property of
the **germline gene**. They are orthogonal — a rearrangement can be in frame with no stop codon
while using a pseudogene V. Calling the first one "functional", as this package did, guaranteed the
confusion. Note the word is already overloaded a third way inside this very package:
`vdjtools.model` trains its recombination models on `LABELS = ("functional", "nonfunctional")`
reads, where *nonfunctional* is exactly the population `filter_productive` removes — and it wants
them, because a rearrangement that never met selection is the cleanest read of recombination.

`filter_productive` **reads the file's own annotation in preference to re-deriving it**: the AIRR
`productive` column if present, else `stop_codon` + `vj_in_frame`, else the `junction_aa` string.
`productive_mask()` returns the predicate *and* which evidence it used, because a caller silently
disagreeing with their own annotation is the failure this prevents. The derived fallback cannot see
a defect in a splicing site or a regulatory element, which AIRR's `productive` can.

`filter_functional_genes` uses IMGT functionality that was already in `model/reference.py` and had
never been wired to a filter. An unresolvable gene call is **kept**: an unrecognised name is a
vocabulary gap, and dropping those rows would report a nomenclature bug as biology.

`filter_functional(keep="coding")` still works and warns.

### Added — `read(..., recompute_frequencies=False)`: use the frequencies as in the file

This was **not possible before**. `_AIRR_ALIASES` had no `frequency` entry — the comment said
outright *"frequency is always recomputed from counts, never read from source"* — and all three
readers passed `recompute_freq=True` unconditionally (`io/read.py`:107, :179, :278). So a
UMI-corrected frequency, or one normalised against anything other than the row counts, was
silently replaced by `count/total` at read time and could not be recovered at any layer above.

`frequency` now aliases from the source, and `read`/`read_vdjtools`/`read_airr`/`read_parquet`
take `recompute_frequencies` (default `True`, unchanged behaviour). A file with no frequency column
derives one either way — there is nothing to preserve. Where `read_airr` collapses rows to
clonotypes a preserved frequency is **summed**, because two rows becoming one clonotype contribute
additively to its share; the legacy converters build `frequency` themselves and are unaffected.

### Added — `filter_length`, and `recompute_frequencies` on every filter that takes it

`filter_length(df, min_len=5, max_len=60)`, both bounds **inclusive**. A data-sanity bound, not a
biological claim: below 5 aa a junction cannot span the Cys104..Phe118 anchors with any diversity
between them, above 60 aa it is beyond what the germline produces. Note these are `junction_aa`
lengths, two residues longer than the IMGT CDR3.

**Frequencies as they are in the file are what you get.** Reading never renormalises; a frequency
is only recomputed when a filter removes rows, and `filter_productive`, `filter_length` and
`filter_functional_genes` all expose that as `recompute_frequencies` (default `True`).

Documented honestly: the switch is undone by the next filter in a chain, because `filter_frequency`,
`filter_segment`, `filter_by_sample` and `downsample` renormalise unconditionally and `select_top`
spells it `renormalize`. Do the frequency-preserving filter last.

### Added — a `Data pre-processing` documentation page

Reading and format conversion, the three axes, length and frequency and segment filtering, error
correction, depth normalisation, pooling and batch correction, and the CLI for all of it. Its
absence is why the axes were conflated in the first place.

### Changed — an unparseable junction is an error, not something to filter away

`sanitise()` used to drop everything that was not a plain amino-acid string. That conflated two
different things, and only one of them is a filtering decision:

- **non-functional** — a stop codon `*`, or the legacy out-of-frame marker `_`. A real
  rearrangement that encodes no receptor. Dropped, and its weight is what
  `vsig:qc:*:nonstd_aa_frac` reports.
- **corrupt** — an ambiguity code (`X`/`B`/`Z`), a lowercase residue, an empty string, a stray
  character. Not a category of receptor: the table is damaged, or was written by something that did
  not agree on the alphabet.

Silently dropping the second hid a broken input *and* inflated the reported non-functional fraction
with junk. `sanitise(df, *, strict=True)` now raises on it; `strict=False` restores the old
behaviour for a corpus known to carry ambiguity codes. A `null` junction stays a drop — absence is
not a damaged character, and a partly-annotated table should still yield a vector.

The character class is not a guess. Measured across **6,047,716 rows** of one bulk-RNA-seq AIRR
cohort, all seven loci:

| | rows |
|---|--:|
| plain 20-amino-acid junctions | 5,600,475 |
| containing a stop codon `*` | 447,241 |
| containing `_` | **0** |
| anything else | **0** |

So the strict default costs nothing on well-formed data and is only ever reached by a file with a
real problem.

### Added — `vsig(..., prefiltered=True)`, and a warning when it should have been passed

A collaborator who removes non-functional rearrangements upstream got a vector that differed from
ours, and the natural reading — a total-frequency denominator computed on the unfiltered table —
is wrong. Every block is computed on the rows that survive `sanitise`, and `work_frame` overwrites
`frequency` with `log2(1+count)/Σ` over those survivors without ever reading the input file's
`frequency`, so an upstream renormalisation cannot propagate.

What differs is that `sanitise` also *reports the weight fraction it dropped*, via
`logit(nonstd_frac, raw.height)`. Pre-filtering drives the numerator to zero **and** shrinks the
denominator, so it moves twice. Measured on 1,168 bulk-RNA-seq blood samples at `tier="standard"`:

| | columns |
|---|--:|
| compared | 688 |
| moved | **7** — one `nonstd_aa_frac` per locus |
| bit-identical | 681, including all 528 geometry columns |
| moved under `compact` / `transfer` / `classify` | **0** |

That column emitted `logit(0, n) ≈ −13.9` on a pre-filtered input: a confident number meaning "this
repertoire is exceptionally clean" where the truth is "we cannot tell". It violated this package's
own rule that holes are load-bearing. `prefiltered=True` now reports `nan`.

`vsig` also **warns** when every present locus reports exactly zero non-functional weight. A single
locus reaching zero is ordinary — IGK does it in 84 of those 1,168 samples, 54 of them with 100+
clonotypes — but every locus at once happens in **0 of 1,168**, the minimum being 5 of 7. So the
detector is on the conjunction, which has no false positives in that corpus, and `prefiltered` is a
parameter rather than something inferred.

### Why this is worth a release rather than a footnote

Costed on real endpoints, the filtering choice moves a **stage-IV immune-toxicity AUC by 0.073** and an
OS C-index by 0.017 for arms that carry the `qc` block. For arms under any `recommended` preset it
moves nothing at all, on 30 of 30 (endpoint, arm) cells.

## 3.9.3 — 2026-08-16

### Changed — the release gate no longer runs the ten slowest tests

The publish gate added **~25 minutes of wall time to every release**, and the comment introducing
it in 3.9.2 asserted the opposite: *"running in parallel with the wheel builds, costs no extra wall
time."* Measured on the 3.9.2 release run (33m03s end to end), that was wrong by a wide margin:

| job | duration |
|---|--:|
| build-sdist | 0m16s |
| build-wheels macos | 2m42s |
| build-wheels ubuntu | 3m07s |
| build-wheels windows | 6m06s |
| **test (the gate)** | **31m18s** |
| publish | 0m43s |

Wheels finished at +6m06s; `publish` started at +32m19s. For comparison 3.9.1, before the gate
existed, took 15m26s. Inside the gate the split is **51 s** of install-and-compile against
**30m21s** of pytest — roughly 1:180, so the C++ build was never the problem.

The obvious fix does not work: `slow` is already deselected by `addopts`, so those 19 tests **never
run in CI at all** and the 30 minutes is entirely non-`slow` tests. Hence a second marker.

**`heavy`** — ten tests, measured at **318 s of a 322 s** three-file run:

| test | measured | why |
|---|--:|---|
| `test_gene_prior.py` (all 4) | 162.9 s | four full `infer_native` EM runs over the 89-allele bundled TRB locus; the C++ E-step is thread-parallel, so it degrades further on a 4-vCPU runner |
| `test_viterbi.py` `[TRB]` (4) + chosen-D | 130.2 s | the pure-Python reference scenario DP, documented at ~600× slower than native; the `[TRA]` twins cost 0.4 s because VJ has no D enumeration |
| `test_dynamics_paired.py::test_pvalues_are_calibrated…` | 25.4 s | a 200k/2M-read replicate pair plus a mandatory negative control |

`publish.yml` now runs `-m "not slow and not heavy"`. **Nothing loses coverage:** `ci.yml` still
runs the full suite on every push and PR at the same SHA — only the release path is trimmed.

Unlike `slow`, `heavy` is **not** in `addopts`, so a plain `pytest` locally still runs all ten.

Measured locally (M3, `HF_HUB_OFFLINE=1`):

| suite | result | wall |
|---|---|--:|
| full (`-m 'not slow'`, what the gate ran) | 1150 passed, 12 skipped, 19 deselected | 387.23 s |
| release path (`-m 'not slow and not heavy'`) | **1140 passed, 12 skipped, 29 deselected** | **72.65 s** |

**5.3× faster**, ten more tests deselected, same pass count otherwise. CI runs ≈4.7× slower than
this machine, so the projected release-gate job was **~6 min against the measured 31m18s** — at which
point it really would finish alongside the 6m06s Windows wheel build, as 3.9.2 wrongly claimed it
already did.

**Measured on this release run, and it beat the projection:**

| release | gate job | wall |
|---|---|--:|
| 3.9.2 | `Test before publishing` | 31m18s |
| 3.9.3 | `Test before publishing` | **4m04s** |

**7.7× on CI**, against 5.3× locally and the ~6 min projected. The runner gains more than the
arithmetic predicted because `test_gene_prior.py`'s thread-parallel C++ E-step degrades worst on a
4-vCPU runner, and that is exactly what `heavy` removes. The gate now finishes well inside the
6m06s Windows wheel build, so it genuinely costs no extra wall time — which is what 3.9.2's comment
asserted without measuring.

Coverage is unchanged: `ci.yml` ran the **full** suite green on this same SHA across ubuntu and
macOS × Python 3.10/3.12 (run 31970256603).


### Added — `Manifest.builder_version`

Which vdjtools built a model, set by `data.build_model` and persisted in `manifest.json`. A germline
defect lives in the **builder**, not the schema, so `model_version` could not answer *was this built
before or after the fix* — for the 3.9.1 J-anchor defect the answer had to be reconstructed by
comparing a shipped model's posterior against an unaffected reference fit. Now it is a field lookup.

Backward compatible: a manifest written earlier reads `""`, which means *predates 3.9.2*, not
*missing*. The bundled models are unchanged in this release and so all read `""`. A new
`test_collapse.py` check asserts the seven `learned` loci always report the **same** builder — a
*partial* regeneration is the dangerous state, because half the set would silently answer a
different question from the other half.

### Fixed — the PyPI upload was not gated on a green test suite

3.9.1 uploaded to PyPI while CI was still running its Test step. Same commit, so the code was
covered, but nothing enforced the **order**, and a CI-only failure could have landed on an
already-published version. `needs:` cannot reference another workflow, so `publish.yml` now carries
its own `test` job — same extras as `ci.yml`, so the OLGA oracle suite actually runs instead of
silently skipping — and `publish` needs it. It runs alongside the wheel builds, so it costs no
extra wall time.

### Note — the bundled models are NOT regenerated, deliberately

3.9.1 fixed `from_olga(derive_orf=True)`, which had reconstructed J CDR3-region germlines as
`full[anchor:]` (the **V** convention) instead of `full[:anchor + 3]`. Fixing the builder does not
fix already-built parquet, so the bundled `learned` human TRA model still carries 11 J germlines
from the wrong side of the anchor, and `TRAJ35` — the one functional allele among them — reads
**~774× low** against the unaffected `arda` fit.

That is left in place on purpose. **0** of VDJdb's 30,937 human TRA records (27,272 unique
junctions, 49 TRAJ genes) use any of the 11 affected alleles — verified as real absence, not a
name-matching artefact, since neighbouring `TRAJ34` (891 records) and `TRAJ36` (494) are well
covered. Regenerating the set is a multi-hour all-loci job, and there is no measured consumer of
the difference. The affected alleles stay pinned by name in `test_collapse.py`.

**If you use human TRA and care about `TRAJ35` usage, use `load_bundled("TRA", "arda")`** — that set
is built from arda germline, where `derive_orf` never runs, so it cannot carry this defect.

Measurements: `bench/results/vdjtools_germline_pgen_shift.md` in the benchmark repo, which also
records what the **3.9.1** collapse fix already moved — 6,676 of 41,322 VDJdb human TRB junctions
(16.2%) went from `Pgen` exactly `0.0` to positive, all of them `TRBJ2-7`, with max
`|Δlog10| = 0.000000` across every junction that was already non-zero.

## 3.9.2 — 2026-08-16

### Added — `Manifest.builder_version`

A germline defect lives in the **builder**, not the schema, so `model_version` could not answer
"was this built before or after the fix" — the 3.9.1 J-anchor defect had to be sized by comparing a
shipped model's posterior against an unaffected reference fit. Set by `data.build_model`, and
backward compatible: `""` means "predates 3.9.2", not "missing".

`test_the_learned_set_is_never_a_mix_of_builders` asserts the seven loci agree on one builder
version. It deliberately does **not** require every model to be stamped — that could only be
satisfied by regenerating, which would encode "always rebuild" as policy. The mixed state is the
dangerous one, and this caught exactly that: an interrupted all-loci run left 4 loci at `3.9.2` and
3 at `""`.

### Fixed — the PyPI upload is gated on a green suite

3.9.1 uploaded to PyPI while `ci.yml` was still in its Test step. Same commit, so the code was
covered, but nothing enforced the **order**, and a CI-only failure could have landed on an
already-published version. `needs:` cannot reference another workflow, so `publish.yml` gets its own
`test` job (ubuntu, 3.12, the same extras `ci.yml` uses, so the OLGA oracle suite actually runs) and
`publish` needs it. It runs in parallel with the wheel builds, so it costs no extra wall time.

### Changed — the bundled models are deliberately NOT regenerated

3.9.1 fixed the builder, and fixing a builder does not fix already-built parquet, so `learned` human
TRA still carries 11 J germlines from the wrong side of the anchor and reads `TRAJ35` ~774x low.
Nothing consumes that: **0 of VDJdb's 30,937** human TRA records use any of the 11 alleles, which is
real absence rather than coincidence — neighbours `TRAJ34` (891) and `TRAJ36` (494) are well covered.
Regeneration is a multi-hour all-loci Aldan-3 job and is not run to populate a metadata field. Users
needing `TRAJ35` usage take `load_bundled("TRA", "arda")`, which cannot carry the defect.

The rule this produced is in `CLAUDE.md` as a decision test rather than a judgement call: **if no
germline entering a locus changed, that locus's model cannot change, so regenerating it is pure
cost.** Measured, not assumed — rebuilding on the same corpus reproduced TRB and TRG
**bit-identically** (max `|dlog10| = 0.0000`), because neither had an affected allele. Only TRA
moved.

## 3.9.1 — 2026-08-16

### Fixed — a collapsed gene could be represented by a **non-functional** allele, making Pgen a silent zero

`collapse_alleles` picks one allele per gene, keeps its germline and relabels it `gene*01`. It
ranked candidates by CDR3-region germline length, then usage, then the allele **name** — but usage
was only ever passed for **V**. For J and D the key therefore degenerated to the name after the
length tie, and `max` took the lexicographically last allele, which was then relabelled `*01`.

Human `TRBJ2-7` is the worst case. `*01` is functional and templates `SYEQYF`; `*02` is an IMGT
**ORF** and templates `SYEQYV`. Both are 19 nt, so length did not separate them and `*02` won on
the name — the collapsed model shipped `*02`'s germline under `*01`'s label:

```
Pgen("CASSIRSSYEQYF" | TRBV19*01, TRBJ2-7*01), bundled human TRB `olga` model
             collapse=True   collapse=False
  3.9.0        0.0             4.7889e-08
  3.9.1        4.9986e-08      4.7889e-08
```

**Exactly zero, with no error raised.** Measured on 864 `TRBJ2-7` nucleotide junctions drawn from
the model itself (`generate(seed=7)`, n=4000): **864 of 864 scored `pgen_nt == 0`** against the
3.9.0 collapsed model, 3 of 864 against 3.9.1 — and those 3 are genuine `*02`-germline draws, which
a single-germline collapsed model cannot represent by construction. Summed over the 864: `0.0` at
3.9.0, `3.357e-07` at 3.9.1, against `3.725e-07` uncollapsed (mean `|Δlog10|` 0.157). This affects
**every released version that shipped these bundled models**, on the default `collapse=True` path.
`load_bundled(..., collapse=False)` reproduces the old numbers for anyone who needs them.

The representative is now ranked by **length → IMGT functionality (`F` > `ORF` > `P`, read from
arda's `cdr3_anchors.tsv`) → usage → prefer `*01` → name**, and the J and D usage marginals are
actually passed. Length still leads: over all 23 bundled models the two orders differ on exactly one
gene, human `TRBV23/OR9-2`, where the only non-pseudogene allele has an **empty** CDR3 germline and
leading with functionality would install it — trading one silent zero for another.

**65 genes change their collapsed germline** (every change length-preserving, 63 of them onto
`*01`): 10 in `olga`/human, 9 in `learned`/human, 29 in `arda`/human, 17 in `arda`/mouse. Beyond
`TRBJ2-7`: `IGKJ2`, `IGKJ4`, `TRAJ47` in all three sources; `TRAJ24/32/37/41`, `IGHJ4`, `IGHJ6`,
`IGKV1-39`, `IGKV3D-20`, `IGHV2-70`, `IGLV1-41` and 7 `IGHD` genes in `arda`; 12 mouse `TRAV` plus
`TRBJ1-1`/`TRBJ1-5`. **Any Pgen, generation or scoring output on a collapsed model can move.**

### Fixed — `from_olga(derive_orf=True)` rebuilt J germlines from the wrong side of the anchor

The CDR3 region lies on opposite sides of the conserved-codon anchor per segment: V runs
`full[anchor:]`, J runs `full[:anchor + 3]`. `_genomic_table` used the V slice for both, so every
ORF/P J allele it reconstructed got the framework *downstream* of Phe118 — still a plausible
in-locus sequence, hence silent. 11 alleles in the bundled `learned` human TRA model carry it
(`TRAJ1*01`, `TRAJ2*01`, `TRAJ19*01`, `TRAJ25*01`, `TRAJ35*01`, `TRAJ51*01`, `TRAJ55*01`,
`TRAJ58*01`, `TRAJ59*01`, `TRAJ60*01`, `TRAJ61*01`); the builder is fixed, but clearing the shipped
parquet needs a model rebuild, so the new test pins those 11 as a named exception.

### Added — permanent anchor tests over every shipped model

`tests/python/test_collapse.py` now asserts, across all 23 bundled models (`olga`/`learned`/`arda`
× human, `arda` × mouse TRA/TRB) and both `collapse=True` and `False`, against arda's per-allele
`functionality`/`status`/`templated_aa` rather than a guessed convention:

- every functional J germline translates, in its anchor frame, to a terminal **F or W**. A
  collapsed row is held to its *gene's* standard, not to whichever allele supplied the germline —
  judging it by that allele would have exempted `TRBJ2-7*02` as "an ORF" and waved the defect
  through. Named exceptions: `TRAJ35*01` (arda records `templated_aa=IGFGNVLHC` at `status=ok`;
  IMGT still calls it F) and `IGHJ6*02` in OLGA's namespace only (OLGA ships it 1 nt short of the
  Trp118 codon; arda's own copy templates `YYYYYGMDVW` and passes);
- no gene is represented by an ORF/pseudogene allele where a functional one exists — 18 genes
  failed this before the fix;
- the CDR3-region germline sits on the documented side of the anchor;
- regression pin: `Pgen("CASSIRSSYEQYF" | TRBJ2-7*01) > 0` in all three sources, with collapsed and
  uncollapsed agreeing to within 0.3 in `log10`.

## 3.9.0 — 2026-08-16

### Added — `Pgen` of a degenerate motif

`native.pgen_aa_degenerate(model, allowed, v=None, j=None)` and `pgen_aa_degenerate_batch` expose
the masked transfer-matrix DP that already backed `pgen_aa` and `pgen_aa_hamming1`. `allowed` is one
entry per position, each a string of permitted residues; `""` or `"X"` means any residue.

This makes the total generation probability of a V/J/length-pinned motif — a VDJdb cluster PWM, say
— a single exact call, with no enumeration and no inclusion–exclusion. Motivated by epitope
precursor-frequency estimation (`appendix/pgen_motif.md`).

Note `pgen_aa` itself still scores an `X` as **0.0**, because `mask_for_aa` matches the genetic code
by exact character. Use `pgen_aa_degenerate` when a position is meant to be a wildcard.

## 3.8.0 — 2026-08-15

Single-cell interop: vdjtools now sits inside the downstream single-cell ecosystem instead of
ending at its own frame.

### Fixed — `paired_pgen` returned nothing but nulls on real CellRanger data

The bug that mattered most here, and it was silent. CellRanger reports **gene**-level V/J calls
(`TRBV10-3`); the model is keyed by **allele**. `native.pgen_aa` raises on a gene name on
purpose — the old `-1` fallback meant *marginalise over every allele* and once returned a Pgen
**2.38x too high** with no error — but `sc.pgen._chain_pgen` caught that with a bare
`except Exception: return None`. Net effect: on the single most common real input, `pgen_alpha`,
`pgen_beta` and `pgen_paired` were **100% null**, with nothing to indicate why. Measured on the
public dCODE donor-4 run: **27,268 of 27,268 receptors null**.

`paired_pgen` now resolves a gene to its representative allele (`*01` where the model has it)
before scoring — deliberately and documented, *not* by falling back to marginalising. Same
dataset: **24,325 of 27,268 now scored** (median paired Pgen 2.1e-19). `resolve_genes=False`
restores exact-allele-only matching, and an all-null locus now emits a `UserWarning` instead of
shipping a silent column. The `except` is narrowed to `(KeyError, ValueError)` with a
non-`str` junction guard, so unrelated failures stop being swallowed.

### Fixed — a barcoded AIRR table was silently collapsed into a bulk repertoire

`io.sniff_format` had no `cell_id` branch, so CellRanger's `airr_rearrangement.tsv` (or any
barcoded AIRR table) sniffed as `"airr"`/`"arda"` and `read_airr` pooled reads **across cells**,
dropping the barcode with no error. It now sniffs as `"airr_cell"` and `io.read` refuses it,
naming `sc.read_airr_cell` instead; `fmt="airr"` still pools on purpose.

### Added — one interchange format, four ecosystems

scirpy, dandelion and scRepertoire all read the same thing: a flat AIRR Rearrangement table with
`sequence_id` + `cell_id`. So `vdjtools/sc/airr.py` is one emitter (`to_airr`) and one inverse
(`from_airr`), and each bridge is a thin adapter — `write_airr`, plus `write_screpertoire`
(`format="airr"|"10x"`). It reconciles the two spellings that otherwise bite: AIRR says
`junction` where vdjtools says `junction_nt`, and scRepertoire's parser reads `consensus_count`
where scirpy and dandelion prefer `umi_count`, so both are emitted.

- **scirpy / scverse** — `to_scirpy` (scirpy's `obsm["airr"]` awkward layout, `index_chains` run
  by default; `gex=` gives a `MuData`) and `from_scirpy`. Writing **delegates** to
  `scirpy.io.read_airr` so no copy of their schema can drift here; reading is ours and needs only
  `awkward`, so consuming someone else's AnnData costs no scirpy install.
- **dandelion** — `to_dandelion` / `from_dandelion`, plus `read_h5ddl`: `.h5ddl` is plain HDF5, so
  a dandelion result opens with `h5py` alone.
- **`push_obs`** — attach vdjtools-computed columns (`pgen_paired`, mispairing flags) to an
  `AnnData.obs` or `Dandelion.metadata` you did not build. Refuses a multi-pair frame rather than
  silently picking one row per cell.

### Added — ingestion

`read_10x` now accepts `filtered_contig_annotations.csv` as well as `all_contig_annotations.csv`
(one CellRanger writer, one layout) and tolerates version drift — `fwr*`/`cdr1`/`cdr2` are CR6+,
`exact_subclonotype_id` CR4+, `sample` only under `cellranger multi`, and `raw_consensus_id` is
used when present rather than required. `read_arda_cells` reads `arda cells` output
(`.contigs.airr.tsv` + `.chains.tsv`), surfacing arda's own per-chain verdict as `arda_status`
**without acting on it** — arda's call and `resolve_chains`' call are independent answers to the
same question. `productive` joins `SC_COLUMNS` so the emitted AIRR table is schema-valid.

### Added — CLI, docs, example

`vdjtools sc` — `convert`, `pair`, `qc`, `pgen`, and
`export --to airr|scirpy|dandelion|screpertoire|screpertoire-10x|airr-cell`, each exposing the
matching library options: `--fmt`, `--require-cell`, `--require-high-conf`, `--consensus`,
`--locus-pair`, `--resolve`, `--flag-mispairing`, `--max-slaves-per-master`, `--drop-mispaired`,
`--source`, `--condition-vj`, `--resolve-genes`, `--alpha-locus`, `--beta-locus`,
`--index-chains`, `--repertoire-id`. (GEX pairing stays library-only --
`to_scirpy(cells, gex=...)` returns the MuData; a CLI flag for it only added a MuData
serialisation step, which failed in CI environments we could not reproduce.) The input format is sniffed from the
**header**, not the filename — a renamed export still works and a bulk table is refused by name
rather than mis-parsed. `sc pgen` reports `scored N/M receptors`, so a naming mismatch is a
number on screen rather than a column of nulls to notice later. A dedicated
`docs/singlecell.rst` (the `usage.rst` section is now a pointer), and
`examples/single_cell_interop.py`, a marimo notebook running the whole path on dCODE donor 4 —
which is what surfaced the Pgen bug above.

`[sc]` gains `awkward` + `mudata`; a new **test-only** `[interop]` extra carries `scirpy` and
`sc-dandelion` (PyPI name; imports as `dandelion`). CI installs it best-effort and reports
whether the round-trip tests ran or skipped, since their dep chains break on new matplotlib.
The format contract itself (`test_sc_airr.py`) has no optional deps and never skips.

## 3.7.3 — 2026-08-15

Housekeeping. The first PyPI release since 3.7.0, so it carries 3.7.1 and 3.7.2 with it.

### Fixed — the iNEXT bootstrap carried a fallback that could never be taken

`stats/inext.py` guarded its `_core` import in a `try/except` and dispatched between the native
bootstrap and the numpy reference at call time. `_core` is a build-time dependency — an install
without it does not exist — so the `except` branch and `_bootstrap_se_dispatch` were dead, and
`inext_batch` raised its own "requires the native _core extension" for a state that cannot occur.
Both are gone; the import is now local to the two functions that need it, which keeps
`import vdjtools.stats` as light as the guard made it. The numpy `_bootstrap_se` stays, unchanged
— it is the reference the tests compare the native kernel against.

### Removed — three dev-notes files the changelog had already absorbed

`NOTES.md`, `ROADMAP.md` and `SUGGESTED_EDITS.md` recorded the phase narrative from before the
changelog existed, and had been drifting from it since. `CHANGELOG.md` is the release-by-release
record; `CLAUDE.md`'s "Open loops" is what is in flight. `CLAUDE.md`, `README.md`,
`docs/index.rst` and the sdist exclude list no longer point at the deleted files.

## 3.7.2 — 2026-08-14

Documentation accuracy. No code change.

### Fixed — the preset-ranking corpus was described as larger than it is

`vdjtools.signature.presets` and `docs/signature.rst` both said the rankings come from "several
hundred study groups, tens of thousands of samples". Counted, the sweep panel is **182 study
groups over 198 accessions, 14,553 samples** — which the same page already stated correctly two
sections earlier ("14,553 samples × 1,369 columns, 182 studies"), so the file disagreed with
itself. Both places now carry the counted figure. The accession list is published in the analysis
repo's `heldout/signature_studies.tsv`, so a reader can check the claim rather than take it.

## 3.7.1 — 2026-08-14

Audit pass. Two fixes, both cases where a feature was reachable from Python and not from the
command line that ships it.

### Fixed — `keep=` stopped at the readers, so the SHM block could never be computed from the CLI

3.7.0 added `keep=` to `read_airr` / `read_vdjtools` / `read_parquet` so `v_identity` — the one
field the signature needs that the canonical eight columns do not carry — could reach
`vsig:shm:IGH:mean_v_identity`. The dispatcher `io.read` and the batch mapper `io.map_samples`
did not take the argument, and those are what the CLI uses: `vdjtools signature` on a file
carrying `v_identity` reported `mask:IGH:shm = 0` and `mean_v_identity = nan`. Both now take
`keep=`, and the `signature` command passes `("v_identity",)`.

The column now populates on any input that has the field. Nothing else changes: `keep=()` is the
default everywhere, and the legacy converters, which narrow to the canonical schema, ignore it.

### Fixed — a coverage level of exactly 1.0 warned its way to the right answer

`mir.signature` passes `cstar = 1.0` deliberately, as an "unreachable" sentinel, for a locus
where no coverage level could be established — the diversity block is then supposed to fail its
own estimability check and mask out. It did, but `_invert_coverage` got there by evaluating
`log(1 - 1.0)` and doing inf arithmetic, emitting three `RuntimeWarning: divide by zero` per
call into the user's terminal. It now returns `inf` directly. Same `m`, same method, same mask —
without the noise.

## 3.7.0 — 2026-08-14

### Added — `vdjtools signature` on the CLI, with the help text as the primary documentation

The command a collaborator actually runs. No Python:

```bash
vdjtools signature --preset classify -m metadata.txt --base-dir samples/ -o sig.tsv
vdjtools signature --preset compact a.tsv b.tsv.gz -o sig.tsv
vdjtools signature --preset classify --describe     # the columns, reading no input
vdjtools presets                                    # the named feature sets, ranked
```

`--help` on both commands carries worked examples, the three `recommended` presets and when each
applies, the pointer to `mir signature` for the geometry half, and the CDR3-vs-junction trap
(a file carrying only IMGT `cdr3_aa` is two residues short everywhere, which shifts the length,
k-mer and Pgen features). `docs/signature.rst` opens with the same quickstart, and the README and
`examples/README.md` lead with `--preset classify` rather than a `specific`-ranked set.

Because that help text is written for a terminal — indented example blocks, which are not valid
reStructuredText — `signature` and `presets` are excluded from the `vdjtools.cli` autodoc, with a
note on the API page saying where to read them instead.

### Fixed — 221 signature tests were invisible to CI

The seven new test files landed in `tests/` while `testpaths = ["tests/python"]`, so a plain
`pytest` collected 789 of 1010 and the CI job (`pytest tests/python -q`) never ran one of them.
Moved into `tests/python/` with the rest; default collection is 1004 passed / 6 skipped.

### Added — `vdjtools.signature`: VSIG, the statistics half of a portable repertoire signature

One repertoire in, a **fixed, named, positional** feature vector out — the object you hand a
collaborator so their matrix and yours are the same coordinate system. The geometry half lives in
`mir.signature`; the shared column contract lives here, because mirpy depends on vdjtools and not
the reverse, and two copies of a contract are not a contract.

```python
from vdjtools.signature import vsig, vsig_cohort, columns, describe
v = vsig({"TRB": df}, tier="standard")
describe("standard")            # column, sig, block, locus, feature, tier, transform, flags
```

Four modules. `layout` is the contract — loci, the `core ⊂ standard ⊂ full` tiers as exact
**index subsets** of one frozen column order, and a per-feature (not per-block) transform
declaration, because a clonality block legitimately mixes a CLR-transformed composition with a
logit-transformed proportion. `transform` is the variance-stabilising layer. `blocks` computes.
`assemble` puts them in order.

Every transform choice is **denominator-aware**, because the alternative silently lies about
shallow samples: Haldane–Anscombe `logit` so `0/3` and `0/500` are different numbers, Anscombe
`arcsine` so a share is defined at exactly zero, and `clr` over the *whole* composition before any
coordinate is selected — shipping *k−1* parts, since all *k* are linearly dependent and would put
a guaranteed zero eigenvalue in any PCA.

Diversity is compared at a **frozen coverage level**, and `estimable()` **refuses** rather than
extrapolates. Real repertoires attain Good–Turing coverage 0.24–0.58, so a textbook `C* = 0.95`
puts every sample into extrapolation, where the same statistic inflates roughly tenfold. A hole a
model can see beats a confident wrong number. For the same reason `clonality` is rebuilt from the
coverage-standardised Hill numbers, `1 − ln(¹D)/ln(⁰D)`: the observed Pielou evenness it replaced
drifted 0.510 over a 667× depth range, against 0.023 for the standardised form.

### Fixed — the CLR zero replacement could consume the composition it was correcting

The textbook multiplicative replacement puts `delta = 0.5/m` on each zero part and scales the rest
by `1 − n_zero·delta`. On a *shallow* composition that is bigger than the composition: three
parts, one observed, `m = 1` gives two replacements of 0.5 and scales the one real part to exactly
zero, whose log is `-inf` — a value that then propagates through every downstream reduction. Found
while emitting a real 4,000-sample corpus. The replaced mass is now capped below half the smallest
*observed* part, which is the only property a replacement needs; the cap is inactive whenever `m`
exceeds the number of parts, i.e. everywhere outside that tail.

### Fixed — `pgen_block` reloaded the recombination model on every call

Loading and collapsing a bundled model costs 0.4–1.8 s; the Pgen batch that follows costs ~0.15 s.
A corpus emission therefore spent 80–95% of its time re-reading seven files it had already read,
and a seven-locus sample paid it seven times over. Memoised per locus: **1.54 s → 0.01 s** on the
second call.

Also worth knowing when emitting a corpus: `pgen_block`/`vsig` default to `threads=0`, meaning
*all cores*, which inside your own process pool means every worker claims the whole machine. On a
16-core box, 14 workers took the load average to 227. Pass `threads=1` there.

## 3.6.1 — 2026-08-14

Audit pass. No library behaviour changes.

- **`examples/emerson_cmv_hla.py` did not parse on Python 3.10/3.11.** A backslash inside an
  f-string *expression* (`f"…{meta['hla'].str.contains(r'HLA-A\*02').sum()}…"`) is 3.12-only syntax,
  but the package declares `requires-python = ">=3.10"`. The regex is hoisted to a local.
- Unused imports and multi-import lines cleaned out of `examples/` (`ruff --fix`).
- `[tool.ruff.lint]` now ignores `E702`/`E741`/`E731` — the paired-short-statement style and `l`/`O`
  loop scalars are deliberate throughout the examples and OLGA oracle shims. `ruff check .` is
  green, so a real finding is visible again instead of being buried in 30 style hits.
- Repo cleanup: 1.7 GB of regenerable artifacts removed (`examples/.data` notebook caches,
  `docs/_build`, `examples/__marimo__`, tool caches, `.DS_Store`), plus four worktrees whose
  branches were already fully merged into `master` (`feature/biomarker-cooccurrence`,
  `feature/dynamics`, `chore/vdjmatch-pin`, `claude/trusting-torvalds-65e7b7`). `feature/cdr3-viterbi`
  and `signature` still carry unmerged commits and were left in place.

Verified at this commit: `pytest tests/python` 783 passed / 6 skipped; `sphinx-build -W` clean.

## 3.6.0 — 2026-08-12

### Added — `infer_nt`: the nucleotide CDR3 behind an amino-acid one

A VDJdb record carries `(V, J, CDR3aa)` and no nucleotides, so none of the boundary markup a
repertoire analysis wants is there. `infer_nt` reconstructs all of it:

```python
from vdjtools.model import infer_nt
sc = infer_nt(model, "CASSLGQAYEQYF", v="TRBV5-1*01", j="TRBJ2-3*01")
sc.cdr3_nt, sc.v_end, sc.d_call, sc.d_start, sc.d_end, sc.j_start, sc.pgen, sc.margin
```

Two stages. A codon-constrained **max-product DP** over every scenario `pgen_aa` sums — germline
positions pinned to their segment, each free N-region position taking the nucleotide that maximises
`P(nt₁)·∏P(nt_k | nt_{k−1})` under the VD/DJ/VJ dinucleotide model; then a `pgen_nt` re-score of the
survivors, because stage 1 maximises the *joint* `P(nt, scenario)` while the contract is about the
*marginal* `P(nt)`. `pgen` on the result is a real `pgen_nt`.

Stage 1 is **native**: the same Murugan/OLGA `Pi_L·Pi_R` transfer matrix `pgen_aa` already uses,
with `max` in place of the sums and the winning `(V, delV)` / `(J, delJ)` carried through the state
(`native.best_aa_scenarios`). It returns *scenarios*, not nucleotide paths — once the scenario is
fixed, recovering the nt string is one cheap DP over the two insertion blocks, so the sweep carries
no back-pointers. A first cut enumerated the scenarios in Python instead and cost 1.7 s per TRB
CDR3; that implementation is kept as the reference the native one is tested against.

Per-sequence cost, 25 generated productive draws per locus. Leaving the calls out is nearly free,
because the DP sweeps V and J either way:

| locus | ms/seq, V/J known | ms/seq, V/J free |
|---|---|---|
| IGK | **0.41** | **0.48** |
| IGL | 0.46 | 0.50 |
| TRG | 0.78 | 0.68 |
| TRA | 0.86 | 2.43 |
| TRB | 3.02 | 3.20 |
| TRD | 7.07 | 7.78 |
| IGH | 87.14 | 85.88 |

So all 80k VDJdb records — TRA and TRB — take **about 2-4 minutes**. IGH is the outlier at 87 ms:
it carries 60+ D alleles and the D placement loop scales with that library, so plan minutes per
10k there rather than per 80k.

**Measured against the exponential oracle** (generated productive draws, V/J pinned, 10 residues):

| variant | TRG | TRA | ms/seq |
|---|---|---|---|
| one scenario, best codon per residue | 9/25 | 4/19 | 0.06 |
| stage 1 only (`keep=1, n_best=1`) | 21/25 | 15/19 | 0.3 |
| stage 1 + marginal re-score (defaults) | **25/25** | **19/19** | 1.2 |

The cheap shortcut fails because a trim chosen before the codons pins a codon the true optimum
would have trimmed away — the same unsoundness as pinning the germline flanks.

**Three call-input modes**, because annotation tables have all three: one allele (the normal mode);
several, as the comma-separated string an ambiguous AIRR `v_call` carries or as a list; or nothing
at all, where the DP marginalizes over every gene — 0.26 ms against 0.23 ms with both pinned, so
unknown calls cost essentially nothing.

Tandem-D is not enumerated in stage 1 (a single D trimmed to zero length already reaches every
middle, so D-D can only reorder candidates, not add them); the stage-2 `pgen_nt` counts it in full.

### Fixed — EM could relearn a genomically impossible D–J pair

A D can only recombine with a J lying 3′ of it, and in TRB the clusters interleave
(TRBD1·TRBJ1·TRBD2·TRBJ2), so `P(TRBD2 | TRBJ1-*)` must be zero. The learned TRB model had it at
**0.0909**, and OLGA's own TRB gives `P(TRBD2*01 | TRBJ1-6*01) = 0.333`.

The constraint is now applied in **both M-steps before normalization** — a post-hoc patch would be
undone by the very next iteration. `reference.forbidden_dj_pairs` derives the forbidden set from
IMGT cluster numbering, `infer.enforce_dj_order` repairs an existing model, and `check_model` gains
an `impossible_dj_pair` check (`warn` for faithful OLGA imports, which must stay byte-exact for the
Pgen invariant; `error` otherwise).

The bundled `learned` TRB is rebuilt with the constraint: **9 iterations instead of 11**, final
log-likelihood **−33.7575** against −33.7604 — the constraint *improves* the fit.

### Changed

- No emoji anywhere in the repository; the `⛔`/`⚠` markers are now `WARNING:` / `NOTE:`.

## 3.5.0 — 2026-08-12

### Fixed — `generate(seed=)` was not reproducible across processes

`collapse_alleles` used unordered polars `group_by`, so a collapsed model's table row order varied
between processes and the *same* seed drew a different allele. `collapse=True` is the default load
path, so this affected ordinary use. The collapse now maintains order throughout. Expectations
recorded from `generate()` before this release are stale; the bundled models are unaffected, since
they ship uncollapsed and are collapsed at load.


### Added — `vdjtools.model.viterbi`: the argmax side of the Pgen DP

`pgen_nt` **sums** over every recombination that could produce a nucleotide CDR3.
**`best_scenario(model, cdr3_nt, v=, j=)`** takes the **maximum** over the same loops and returns
the single most likely one — which *is* the V/D/J boundary markup:

```python
from vdjtools.model import best_scenario
sc = best_scenario(model, cdr3_nt, v=v_call, j=j_call)
sc.v_end, sc.d_call, sc.d_start, sc.d_end, sc.j_start   # 0-based, half-open, CDR3-nt space
```

It re-derives nothing: every probability comes from `prepare()`'s tables, and the D placement is a
max-product mirror of `pgen._d_middle` over the same `P(D|J)·P(delD|D)·Pins(VD)·Pins(DJ)` terms.

WARNING: **The D therefore obeys `P(D|J)`.** TRBD2 lies 3′ of the whole TRBJ1 cluster, so deletional
joining can never produce a TRBD2–TRBJ1 pair and the model encodes that as a zero. An earlier draft
here chose D by longest exact substring and ignored `j` entirely — it would have called the
impossible pair. There is a regression test.

Validated on 200–500 generated draws per locus: the V span **is** the V germline, the J span the J
germline, the D span the D germline, `scenario_p ≤ pgen_nt` (a maximum cannot exceed the sum it is
taken over), and `scenario_p` recomputes exactly from the reported path's own table entries.

`infer_nt_bruteforce` is an exact but exponential **oracle for tests**.

NOTE: Both functions assume an **in-frame** CDR3 (`len(nt) == 3 × len(aa)`). Real productive receptors
satisfy this; an out-of-frame draw does not (measured: 8 aa against 25 nt).

### Changed

- `arda-mapper` pinned to **>= 2.19.0** (was >= 2.5.5).

### Fixed — `generate(model, n, seed=)` was not reproducible across processes

Same seed, same wheel, a **new interpreter → a different draw**. Within one process it looked
perfect, which is why no existing test caught it.

Root cause: **`collapse_alleles`** — which `load_bundled(..., collapse=True)` runs by default —
built its tables with unordered polars `group_by().agg()`. `group_by` is a multithreaded hash
aggregation, so the collapsed table's **row order varied per process**; `_cum` then assigned the
same cumulative interval to a different allele, and the same `rng.random()` drew a different one.
`_pick` and `default_rng` were correct throughout — the ordering beneath them was not.

NOTE: **Not hash randomisation.** `PYTHONHASHSEED=0` did not help, which is what ruled it out and
pointed at the aggregation. Same class as the nondeterminism recorded against arda's `correct`
stage.

Fixed by `maintain_order=True` on all 12 `group_by` calls in `collapse.py` and all 7 in
`generate.py`. Verified identical across 5 separate processes on TRA, TRB and IGH.

WARNING: **This changes generated output** for a given seed — it has to, since the old order was
arbitrary. Any recorded expectation from `generate()` predating 3.3.0 must be re-derived.

New tests run the sampler in a **subprocess**, because an in-process test agrees even with the bug
present; they fail 4 of 5 without the fix.

## 3.4.0 — 2026-08-12

Follows 3.3.0's model workshop with the defects that workshop then found, and retrains every
bundled model.

### Added

- **Live EM progress, checkpointing and exact resume.** `infer`/`infer_native` take
  `progress=callable(iter, loglik, rel_change, n_scoreable)` — `infer.print_progress()` is a
  ready-made one — so a long fit is visibly converging rather than merely running; the relative
  change it reports is exactly what is compared against `tol`. They also take
  `checkpoint=DIR`/`checkpoint_every=N`, saving the model after each iteration (written to a
  sibling directory and swapped in, so a kill mid-write leaves the previous checkpoint loadable),
  and `infer.resume(DIR, seqs)` continues from one. **Resuming is exact**: 3 iterations plus a
  resumed 4 give the same log-likelihood *and* the same tables as an uninterrupted 7, and the
  training log spans every attempt. Exposed as
  `vdjtools model learn -v --checkpoint DIR --resume DIR`; `vdjtools model build -v` additionally
  stops swallowing arda's mapping output. This matters because IGH's EM enumerates ~1,225 D pairs
  per read against TRB's 9 — 12 minutes per iteration on a 112-core node, and more than 110 minutes
  without finishing one on a laptop.

### Changed

- **All seven bundled `learned` models retrained** on the full non-functional read corpus (every
  available read; out-of-frame *and* stop-codon, since both escaped selection). Every one converged
  on tolerance with a step-by-step monotone log-likelihood, and `check_model` reports zero errors:

  | locus | clonotypes | iters | log-likelihood |
  |---|---|---|---|
  | TRA | 34,238 | 10 | −21.52 → −20.01 |
  | TRB | 122,703 | 11 | −37.01 → −33.76 |
  | TRG | 14,305 | 8 | −21.90 → −20.33 |
  | TRD | 10,915 | 7 | −40.94 → −35.79 |
  | IGH | 141,607 | 9 | −57.16 → −49.88 |
  | IGK | 256,347 | 6 | −16.18 → −14.28 |
  | IGL | 23,469 | 12 | −6.63 → −6.00 |

  Each now ships its training log, so a model states what it was fitted on. IGH was trained on a
  112-core cluster node; the rest fit comfortably on a laptop.
- **Ambiguous junction bases are substituted, not dropped.** They previously crashed EM with a bare
  `KeyError` from inside the native encoder. Both training entry points now substitute `A` by
  default and warn with the count (`ambiguous=None` drops instead) — it affects ~0.01% of these
  reads, so dropping cost sample size for nothing. It is a substitution, not a marginalization.

### Fixed

- **`collapse_alleles` could give a gene a germline it could not use — in the default path.**
  `load_bundled(..., collapse=True)` picks one representative allele per gene and averages the
  other alleles' conditionals onto it. The representative was chosen by usage alone, but IMGT ships
  some alleles with a **truncated** CDR3-region germline: human `IGKV3-20*02` is 11 nt against
  `*01`'s 30 and carried the higher learned usage, so it became the gene's germline — relabelled
  `*01`, which was doubly misleading — and 25% of the gene's own averaged deletion distribution
  landed on trims the 11-nt germline cannot reach. The Pgen DP never visits those, so that quarter
  of the probability vanished from every Pgen through IGKV3-20 instead of being redistributed.
  The representative is now chosen by **germline length first, usage second** (alleles of a gene are
  near-identical through the CDR3 region, so a large length gap means an incomplete database entry),
  and the collapsed deletion conditionals are **projected onto the representative's reachable
  support and renormalized**. Every bundled model is now clean at `collapse=True`.
- **A failed `arda` run reported nothing but an exit code.** `annotate_reads` passed
  `capture_output=True` and let `CalledProcessError` propagate, so arda's own message was swallowed
  — precisely how the CLI rename in 3.3.0 stayed invisible. It now raises with arda's stderr and the
  installed `arda-mapper` version.
- `reference.read_fasta` wraps arda's FASTA parser with gzip support; arda's opens with plain
  `open()`, so a `.gz` reached it as mojibake and died on the first byte.

### Documented

- **Known quirks of the OLGA models** (`docs/model.rst`). The bundled `olga` set is a bit-faithful
  import, and that includes its defects. Verified against OLGA's raw `model_marginals.txt` with
  OLGA's own parser: deletion mass on unreachable trims is **OLGA's**, with identical fractions to
  4 dp (`IGHV4-30-4*01` 100% — its Pgen is identically zero in OLGA too — `IGKJ4*02` 80.9%,
  `TRAV20*03` 54.7%), and our Pgen matches olga-pip exactly (ratio 1.000000) on every sequence OLGA
  will score. Correcting it would break the exact-OLGA-Pgen invariant, so `check_model` reports it
  at `warn` rather than `error` for an OLGA-sourced model. Also covers the empty-germline ORF genes,
  protocol-specific V/J usage, and OLGA's refusal of out-of-frame input.

## 3.3.0 — 2026-08-12

The recombination model becomes a workshop: buildable on your own reference, checkable, comparable,
scoreable and extendable. See the new [user guide](https://docs.isalgo.dev/vdjtools/model.html) and
`examples/model_workshop.py`.

### Added

- **Custom V(D)J reference libraries.** `model.io.from_germline(germline_df, locus=...)` builds a
  model on any germline library; `from_arda` is now a six-line wrapper over it (its output is
  byte-identical, verified table by table on TRA/TRB/TRG). `reference.read_germline_fasta(v, j, d,
  anchors=)` reads your own FASTA — segment comes from which argument a file is passed as, so no
  header convention is assumed — and `reference.validate_germline` audits it. The audit includes
  the two anchor-frame checks (V starts on a Cys codon, J ends on Phe/Trp) that catch the most
  damaging custom-library mistake: a CDR3 anchor one codon off shifts every deletion profile by a
  constant and nothing downstream complains.
- **`model.check.check_model`** — a consistency audit returning a tidy issue frame
  (`severity check event segment allele detail value`) rather than raising, so every problem in a
  model is visible at once. Covers normalization and probability range, alleles missing from the
  germline (or absent from the marginals), functional genes stuck at `P = 0`, unreachable deletion
  mass, incomplete dinucleotide tables, and VDJ/VJ event-set mismatches. Deletion reachability is
  derived from the Pgen DP (`ndel = len(cut_segment) − contributed − max_palindrome`) and reported
  as the **fraction** of each allele's mass that is lost, ranked — up to 25% for `IGKV3-20*01` in
  the bundled IGK model. `vdjtools model check` exits 1 on any error-severity issue.
- **`model.score`** — likelihood and diversity. `model_fit` reports log-likelihood, free parameters,
  AIC and BIC; it uses **nucleotide** Pgen, because `Σ Pgen_nt = 1` makes the log-likelihood proper
  (amino-acid Pgen sums only the in-frame, stop-free fiber, so its missing normalizing constant
  differs between models and it is a relative score only). A sequence the model cannot generate is
  counted in `n_scoreable`, never turned into `-inf`. `free_params` counts **occupied cells**, not
  rows, and drops undefined and unreachable conditional groups — the difference between ~700 and
  ~3,600 parameters for human TRB's `v_3_del` alone. Also `pgen_frame`, `compare_pgen` +
  `pgen_summary` (KS, Spearman, and the headline one-sided coverage counts), and `pgen_spectrum`.
- **Information content and total diversity.** `analyze.total_entropy` gives each recombination
  event's contribution to the scenario entropy (the dinucleotide term is `E[length] × H_step`), and
  `score.diversity` adds a Monte-Carlo sequence entropy with its standard error plus both Hill
  numbers — `2^H` and `1/E[Pgen]`. Human TRB: ~52 bits per rearrangement, ~45 bits per sequence,
  ~3·10¹³ effective sequences.
- **Model comparison.** `analyze.compare_models` reports per-event total variation, `tv_max` and
  Jensen-Shannon over the union of both models' realizations with zero fill, weighted by the
  parent's marginal; `by="gene"` bridges different germline namespaces, and an event factorized
  differently in the two models is flagged `schema_differs` rather than joined. Plus
  `compare_usage` and `compare_net_dot` (a bnlearn `compare_networks`-style graph).
- **Training log.** Every EM fit now appends a run to `model.training["runs"]`, persisted beside the
  model as a `training.json` sidecar and readable as a table with `infer.training_frame`. Both the
  field and the sidecar are optional, so every previously-saved model still loads (with
  `training is None`).
- **`infer.infer_frame`** fits from a clonotype frame, building the per-read V/J masks for you, and
  **`infer.extend_alleles`** adds alleles from a larger germline library.
- **`data.build_all`** runs the full corpus pipeline — fetch FASTQ, map with arda, collapse, EM —
  parallel across chains, exposed as `vdjtools model build`. Two arda-mapped clonotype examples
  (human TRA 34,238 and TRB 100,000 out-of-frame clonotypes, 2 MB total) ship in
  `tests/python/fixtures/model_reads/` as gzipped FASTA with the V/J/D calls in the header, and
  load offline with `data.load_prepared` — no network, no arda, no mmseqs2.
- **`vdjtools model` CLI sub-app** with 13 subcommands. A model is named as a directory or as
  `LOCUS[:source[:organism]]`.
- **Table export/import**: `marginals_frame` / `set_marginals`, and `save_model(..., fmt="tsv")`
  with format auto-detection on load, so a hand-edited TSV directory is a first-class model input.

### Changed

- **`pgen_nt` now releases the GIL**, the one Pgen binding that still held it after the Phase-13
  batch work. Threaded nucleotide Pgen is **11.7× faster** on this Mac and bitwise-identical to the
  serial result.
- `extend_alleles` **preserves each pre-existing gene's total usage**. Alleles of one gene are
  alternative versions of the same gene — a diploid carries at most two — so a richer library must
  split a gene's mass more finely, never multiply it. Seeding each new allele at its gene's average
  without this correction moved gene-level V usage on human TRB by up to 6 percentage points,
  silently reweighting every Pgen through those genes.

### Fixed

- **`annotate_reads` was calling an arda CLI that no longer exists.** It shelled out to
  `arda rnaseq map -o …`, but arda 2.19 turned `rnaseq` into the full map→assemble→correct preset
  with no stage positional and `-p/--out-prefix` in place of `-o`, so every real invocation exited
  2 — i.e. the whole model-training pipeline was broken against the installed arda. Stage-1 mapping
  is `arda map`; the pin is now `arda-mapper>=2.19.0`. Worth recording that a `--help` smoke test
  would *not* have caught this: typer short-circuits `--help` before argument parsing, so
  `arda rnaseq map --help` still exits 0.
- **Ambiguous junction bases crashed EM** with a bare `KeyError: 'N'` from inside the native
  encoder. `infer_frame` and `build_model` now substitute `A` by default and warn with the count
  (`ambiguous=None` drops the clonotype instead) — see `infer.sanitize_junctions`.
- `reference.read_fasta` wraps arda's FASTA parser with gzip support; arda's opens with plain
  `open()`, so a `.gz` reached it as mojibake and died on the first byte.
- `tests/python/test_io_hf.py::test_control_native_schema_capped` called `list_repo_files` outside
  the `hf` fixture's guard, so an offline run failed instead of skipping.

## 3.2.0 — 2026-08-09

### Fixed

- **Adaptive/immunoSEQ gene names were wrong for 100 of the 161 tokens** seen in the IMMREP25
  release + the pairSEQ mock cohort (22,058 of 44,000 gene calls), and **every one of those 100
  outputs is a gene name absent from the IMGT human reference**. `_adaptive_to_imgt` normalised
  Adaptive tokens with a global `re.sub(r"0([1-9])", r"\1", …)`, which always re-emits the trailing
  `-01` as an IMGT *subgroup*: `TCRAJ39-01 → "TRAJ39-1"` (no human TRAJ gene has a subgroup),
  `TCRBV09-01 → "TRBV9-1"`, `TCRBD01-01 → "TRBD1-1"`. Slash ties (`TCRBV03-01/03-02`), family-only
  calls (`TCRBV20-X`) and co-locus names (`TCRAV38-02`, IMGT `TRAV38-2/DV8`) were passed through
  verbatim. Any consumer resolving gene names against a germline reference silently lost the rows —
  tcrdist3 dropped 100 % of both cohorts.

  Whether a token's trailing group is a subgroup or an allele is a per-family fact (`TCRAV01-01` =
  `TRAV1-1`, but `TCRAV22-01` = `TRAV22`), so no regex can decide it. `read_immunoseq` now resolves
  V, D **and** J calls through a shipped CDR-validated table, `resources/adaptive_imgt_map.tsv`
  (163 tokens; the choice among candidates is decided by exact-matching germline CDR1+CDR2 —
  provenance in `SOURCES.md`, rationale in `appendix/adaptive_imgt_map.md`). Tokens outside the
  table fall back to the legacy rewrite, so unknown input behaves exactly as before.
  The legacy Groovy `CommonUtil.extractVDJImmunoSeq` has the same defect — a v1 bug inherited by
  v2, not a porting regression.
- **`read_mixcr` accepts every MiXcr count spelling** (`cloneCount`, `readCount`,
  `uniqueTagCountMolecule`) — a v4 `-readCount` export used to raise outright.

## 3.1.2 — 2026-07-30

### Changed

- Bumped the `vdjmatch` floor to `>=0.1.2` — a fresh install of 3.1.1 could resolve `vdjmatch`
  0.1.1, whose `cluster.overlap()` raised `SchemaError` on any query with zero fuzzy matches
  (`a_idx`/`b_idx` defaulted to a `Null` dtype with no hits to infer from, then failed to join
  against the `Int64`-typed lookup frames); this broke `vdjtools.overlap.fuzzy.fuzzy_overlap` and
  CI's `test_overlap_fuzzy.py::test_fuzzy_no_match_empty_and_zero_metrics`. Fixed upstream in
  vdjmatch 0.1.2.

## 3.1.1 — 2026-07-30

### Fixed

- **`biomarker.association(match="fuzzy")` re-ran the full-cohort search on every call**, even
  when a caller tests many phenotype designs (e.g. one per HLA gene) against the same
  cohort/key/candidates/scope — the search depends only on those, never on the design. At
  full-corpus scale (~50k donors) this drove peak memory from 48G to 256-350G per SLURM task on
  diverse BCR light chains, entirely from redundant `collect()`/`.to_list()` work repeated once
  per design instead of once. Added `prepare_fuzzy_features(cohort, key, candidates=, scope=) ->
  FeatureFrame` + `association(..., features=)` so the search is opt-in-cacheable: build it once,
  reuse across every design against the same cohort.
- Documented that `level_col`'s memory cost is **multiplicative, not additive**: `association()`'s
  feature-join duplicates every matched row once per design level (correct behaviour — each level
  needs its own incidence table — but easy to miss from the prior wording).

### Changed

- Bumped the `seqtree` floor to `>=0.6.1` — fixes a corrupted Miyazawa–Jernigan A–N contact
  energy in `structural()` (0.6.0) and names the offending sequence/index in `gapblock_matrix`'s
  alphabet error instead of just the bad symbol (0.6.1).

## 3.1.0 — 2026-07-28

### Added

- **The bundled `arda` model set is reachable.** `_bundled/arda/` has shipped 9 EM-refit models —
  the 7 human loci **plus mouse TRA/TRB** — in every wheel, but no public call could reach them:
  `SOURCES` listed only `("olga", "learned")`, and `load_bundled` keyed
  `_bundled/<source>/<LOCUS>` while the arda directories are `<organism>_<LOCUS>`, so
  `list_bundled()` reported them as absent. (`from_arda` is not a substitute — it returns the
  *placeholder* marginals meant to be refit by `infer_native`, not these refit ones.)

  `load_bundled` gains a keyword-only `organism=` (default `"human"`) and derives the directory key
  per set. This is the only bundled set in the **arda IMGT allele namespace** — the frame that
  arda-annotated pipelines such as `mirpy`'s prototypes and baked germline distances live in — and
  the only bundled set covering a non-human organism.

### Changed

- `load_bundled` now **raises** when a non-human `organism` is asked of the human-only `olga` /
  `learned` sets, instead of silently handing back the human model; `FileNotFoundError` lists the
  keys that *are* available for that set. Every existing caller passes `(locus, source)`
  positionally and is unaffected.

### Repository

- Consolidated the example notebooks: the old `notebooks/` directory was merged into
  **`examples/`**, so every marimo explorer now lives under `examples/` (docs / README / skills
  updated to match). Examples are not shipped in the wheel — this is a repository-layout change only.
- Merged the two aging notebooks into one **`examples/aging.py`**: it now covers the
  cohort-streaming stats (`diversity_cohort` / clone-size / spectratype), the coverage-standardized
  iNEXT diversity + rarefaction, and the pairwise-overlap→MDS divergence — the union of the old
  `aging.py` and `aging_airr_benchmark.py` (both removed, along with the now-unused
  `aging_manifest.json`). Fixes a latent `cdr3_aa`-vs-`junction_aa` key bug in the old benchmark
  notebook (stale since the v2.2.0 junction rename).
- **Every example is now a marimo notebook** — converted the two remaining plain scripts
  (`emerson_cmv_hla.py`, `scale_cohort.py`).
- **No user-specific absolute paths in the examples/tests** — each notebook resolves data from a
  gitignored **`./data_dump/`** directory first (symlink your copies there), then falls back to
  HuggingFace. VDJdb is fetched from the latest `antigenomics/vdjdb-db` release (cached to
  `./data_dump/`) instead of a hardcoded local checkout.

## 3.0.0 — 2026-07-18

### Added — longitudinal clonotype dynamics (`vdjtools.dynamics`)

- **Recapture model** (`dynamics.capture`) — the VDJtrack size-bucket model: clonotypes binned
  singleton / doubleton / tripleton / large, Poisson capture probability `P = 1 − exp(−f·R)`,
  `Beta(captured, missing)` credible intervals, and the group-effect test — a log-linear
  `log(recapture) ~ size + group + log(div_ratio)` plus a per-bucket paired t-test across donors.
  Python port of Pavlova, Zvyagin & Shugay, *Front Immunol* 2024
  ([10.3389/fimmu.2024.1321603](https://doi.org/10.3389/fimmu.2024.1321603)).
- **Metaclonotype-grouped testing** (`dynamics.test_metaclonotypes`) — collapse a 1-Hamming
  (`scope="1,0,0,1"`) or 1-Levenshtein (`"1,1,1,1"`) CDR3 ball into one feature before the paired
  test, for power on convergent expansions.
- **edgeR NB-exact caller** (`dynamics.expansion_test`) — TMM normalization + qCML common
  dispersion + the negative-binomial exact test (as a Beta-Binomial conditional); the paper's §2.5
  complementary per-clone caller.

  (Complements the existing per-clonotype `dynamics.test_pair`, Ayestaran 2024.)

### Added — cohort-streaming summary statistics

- `io.map_samples(fn, items, *, workers=)` — thread-parallel per-sample reduce, `O(workers)`-sample
  peak memory, results in input order.
- `stats.diversity_cohort(cohort)` — the whole cohort's diversity table in one streamed
  count-spectrum pass, bit-exact vs the per-sample path.
- A `by=["sample_id"]` group-prefix on `spectratype` / `vj_spectratype` / `segment_usage` /
  `vj_usage` / `kmer_profile` / `v_kmer_c_profile` / `physchem_profile` — the whole cohort in one
  fused `group_by` over a `scan_cohort` LazyFrame.
- CLI `diversity` / `spectratype` / `segment-usage` / `overlap` gain `--threads N` (parallel over
  samples) and `--cohort DIR` (one streamed pass over a pre-ingested Parquet cohort); the `overlap`
  command now pre-aggregates each sample once.

### Added — CLI & packaging

- New `vdjtools` subcommands: **`convert`** (read any supported format — native / AIRR / Parquet /
  MiXcr / MiGec / MiTCR / immunoSEQ / IMGT / Vidjil / RTCR / TRUST4 / arda — and write the canonical
  table), **`downsample`**, **`filter`** (coding / non-coding / frequency / V-J segment), and
  **`pool`** (flat pool or incidence `--join`).
- Every command's `-o` is now **format-aware**: a `.parquet` / `.pq` path writes Parquet, anything
  else (or stdout) writes TSV.
- Development switched to **uv** — one repo-local `.venv`, no conda. `setup.sh` is rewritten to be
  uv-first (with a `python -m venv` fallback) and **portable across bash and zsh**. `environment.yml`
  is now optional, needed only for MMseqs2 (arda's aligner) + the slow arda round-trip tests.

### Added — notebooks (marimo, `[examples]` extra)

- `examples/vaccination_tracking.py` — clonotype tracking + recapture model across YFV / influenza
  / TBE vaccination time courses.
- `examples/aging.py` — cohort-streaming diversity / clone-size / spectratype across the Britanova
  ageing cohort.
- `examples/ankspond_motif.py` — the ankylosing-spondylitis TRBV9 "AS27" motif (disease vs HLA-B27
  carriage; Komech 2018).

### Fixed

- **C++ CI version drift** — the native `version()` is now single-sourced from `pyproject.toml`
  (parsed by CMake into the `VDJTOOLS_VERSION` compile definition), and both the C++ and Python
  version tests assert *agreement* rather than a hand-copied literal. A release bump can no longer
  redden CI (as the 2.9.0 bump did, leaving `tests/cpp/test_core.cpp` asserting `"2.8.0"`).
