# SOURCES

Provenance of every dataset used or produced by vdjtools. Never guess a source — record it here.

## Germline reference (canonical)

**arda's germline library is the single source of germline truth** — every V/D/J germline
sequence and CDR3 anchor resolves from it by allele name (`vdjtools.model.reference.load_germline`),
so annotation ↔ scenarios ↔ stitching ↔ Pgen share one coordinate frame. arda's anchor
convention is byte-identical to OLGA's (0-based Cys104/[FW]118 offset into the full germline).

| Dataset | Origin | Format | Notes |
|---|---|---|---|
| V/J germline + CDR3 anchors | arda `database/vdj/<org>/cdr3_anchors.tsv` (via `arda.cdr3fix.load_anchors`) | TSV, per-allele | CDR3-region germline + anchor + functionality; **full-length V/J germline not shipped** (build-time only) — a P1c/stitching prerequisite |
| D germline | arda `database/vdj/<org>/d_germlines.fasta` | FASTA `>LOCUS\|allele` | full D germline, no anchor |

**OLGA/IGoR model files contribute recombination probabilities only** (and, for bootstrap
models, their own IMGT-vintage germline is kept for exact-Pgen fidelity — see the P1a note).

## Phase 1 — model engine (bootstrap)

| Dataset | Origin | Format | How to obtain | Provenance |
|---|---|---|---|---|
| OLGA default models (7 human loci + mouse) | `mirpy/mir/resources/olga/default_models/{human_T_alpha,beta,gamma,delta, human_B_heavy,kappa,lambda}/` (reference only — do not copy mirpy code) | `model_params.txt`, `model_marginals.txt`, `V/J_gene_CDR3_anchors.csv` | read from the local mirpy checkout | published OLGA models (Sethna et al.); **derived** generative-model parameters |
| Synthetic out-of-frame training seqs (~100k/chain) | OLGA generation via mirpy tooling; sample at `mirpy/tests/assets/olga_humanTRB_1000.txt.gz` | TSV (V call, J call, CDR3 nt) | generate with OLGA from the models above | **computed** (Monte-Carlo draws); **no tandem D by design** — see plan Phase 1 note |
| Real AIRR reads (model training, all 7 loci) | HF dataset [`isalgo/airr_model_read`](https://huggingface.co/datasets/isalgo/airr_model_read) (owner, cc-by-nc-nd) | raw 5'-RACE **FASTQ** `{group}/{CHAIN}.{functional,nonfunctional}.fq.gz`; groups = `human`, `human_fetal` (TdT-low), `mouse` | `vdjtools.model.data.prepare(group, chain, label)` = fetch (lazy `huggingface_hub`) → `arda rnaseq map` (V/D/**D2**/J + junction + productivity) → `unique_clonotypes` (dedup to `(v,j,junction)`) | **experimental**; out-of-frame reads train the models (EM), functional reads are the selection test set. D-D annotated per-read via `d2_call`. `yields.tsv` records per-bucket read counts |

## Golden fixtures (tests)

| Dataset | Origin | Use |
|---|---|---|
| IGoR model files (human/mouse, all loci) | `IGoR-models/` (local) | model-loader + Pgen oracle fixtures; canonical `model_parms`/`model_marginals`/anchors format |
| OLGA `default_models/` | OLGA pip package | Pgen oracle (`pip install olga`), as seqtree CI does |
| Legacy input-format samples | legacy vdjtools `src/test/resources/samples/*.txt.gz` (on the `legacy-1.x` branch) | format-conversion conformance (MiXcr, MiGec, ImmunoSeq v1/v2, ImgtHighVQuest, Vidjil, RTCR, …) |

## Carried-over resource data (from legacy vdjtools, `legacy-1.x` branch)

| File | Use |
|---|---|
| `resources/profile/aa_property_table.txt` | amino-acid physicochemical properties (Phase 3) |
| `resources/profile/cdr3contact.txt` | CDR3 contact-probability estimate (Phase 3) |
