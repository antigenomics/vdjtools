# ROADMAP

vdjtools 2.0.0 — phased rewrite. Gitflow: `master` (tagged releases) ← `dev` (integration)
← `feature/*` (one per phase). Merge a feature to `dev` when its tests are green; cut alpha
tags off `master` as phases land; `v2.0.0` at Phase 8.

Legacy v1.x is preserved on the `legacy-1.x` branch and under all existing tags.

| # | Branch | Scope | Status |
|---|---|---|---|
| 0 | (master root) | Git surgery, GPLv3 relicense, repo scaffold, CI/docs/publish, AIRR+polars IO layer | in progress |
| 1 | `feature/model-engine` | Native V(D)J model: polars marginals (D-D), arda scenarios, contig stitching, EM from OLGA-synthetic bootstrap, native Pgen/generation (pybind11), validation vs OLGA | planned |
| 2 | `feature/repertoire-stats` | Diversity (Chao1/ChaoE/Shannon/Simpson/d50/Efron–Thisted, exact+resampled+rarefaction+quantile), spectratype, V/J/VJ usage | planned |
| 3 | `feature/cdr-features` | CDR physicochemical profiles (regions × properties), k-mer / V+k-mer summaries | planned |
| 4 | `feature/overlap-tcrnet` | Sample overlap + TCRnet via vdjmatch/seqtree; pairwise-distance matrix, clustering/MDS, tracking | planned |
| 5 | `feature/preprocess` | Downsampling, error-correction, decontaminate, filters, VJ-usage batch-effect correction, pool/join | planned |
| 6 | `feature/biomarker-assoc` | Fisher incidence association vs HLA/condition/chain-pairing; metaclonotype grouping | planned |
| 7 | `feature/singlecell-interop` | AIRR Cell / 10x paired-chain interop, single-cell metadata bridge, paired α/β Pgen | planned |
| 8 | `feature/cli-docs` | typer CLI, full Sphinx docs, example notebooks, **v2.0.0** release | planned |

**Design principles**: AIRR schema + polars everywhere, minimal OO. Python-first — native C++
(pybind11) only for the Pgen DP, generation sampler, and EM E-step inner loop. Delegate overlap/
TCRnet to vdjmatch and annotation/markup to arda rather than reimplementing. **arda's germline
library is the single source of germline truth** — all V/D/J germline + CDR3 anchors resolve from
arda by allele name (`model.reference.load_germline`), so annotation ↔ scenarios ↔ stitching ↔ Pgen
share one coordinate frame.
