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

## Phase 7 — single-cell (dCODE)

| Dataset | Origin | Format | How to obtain | Provenance |
|---|---|---|---|---|
| 10x dCODE dextramer donors (1–4) | HF dataset `isalgo/airr_benchmark`, path `dcode/vdj_v1_hs_aggregated_donor{1..4}_{all_contig_annotations,consensus_annotations,binarized_matrix}.csv.gz` | gzipped CSV (10x CellRanger VDJ all-contig / consensus annotations; CITE-seq binarized dextramer matrix: `barcode` + 50 `*_binder` boolean cols) | `huggingface_hub.hf_hub_download(repo_id="isalgo/airr_benchmark", repo_type="dataset", filename=...)` (fetched at test time into the HF cache; skips cleanly offline) | **experimental** — 10x Genomics "A New Way of Exploring Immunity" dCODE dextramer single-cell TCR + surface-marker dataset; antigen labels are the single-True `*_binder` column per cell |

## Phase 6 — biomarker association (Emerson benchmark + validation)

| Dataset | Origin | Format | How to obtain | Provenance |
|---|---|---|---|---|
| Emerson HIP cohort (786 subjects) | HF dataset [`isalgo/airr_hip`](https://huggingface.co/datasets/isalgo/airr_hip) — redistributed from Adaptive immuneACCESS **Emerson-2017-NatGen** | per-subject VDJtools tables `corr/HIP#####.txt.gz` (`count freq cdr3nt cdr3aa v d j VEnd…`); `metadata.txt` (TAB-sep: `file_name sample_id age race sex cmv hla`) | `examples/emerson_cmv_hla.py` → `huggingface_hub.snapshot_download(repo_type="dataset", allow_patterns=["corr/{sample}.txt.gz"])`; ingest with `io.ingest_cohort(fmt="vdjtools")` | **experimental** TCRβ repertoires + phenotypes (Emerson et al., *Nat Genet* 2017, doi:10.1038/ng.3822). `cmv` ∈ {`+`,`-`,`NA`}; `hla` = 2-digit HLA-A/B only (`HLA-A*02`); ⚠ `race` contains commas — split on TAB. No discovery/validation split column |
| VDJdb (CMV validation target) | local checkout `/Users/mikesh/vcs/code/vdjdb-db/database/vdjdb.slim.txt` (canonical `antigenomics/vdjdb-db`, 2024-06 release; 2-digit HLA matches airr_hip) | TSV, 16 cols: `gene cdr3 species antigen.epitope antigen.gene antigen.species … v.segm j.segm … mhc.a mhc.b mhc.class … vdjdb.score` | read with polars; filter `gene==TRB & species==HomoSapiens & antigen.species~CMV` | **curated** TCR↔epitope database; newer 4-digit dump at `/Users/mikesh/vcs/code/vdjdb-iedb-concordance/vdjdb_dump_2026/vdjdb.slim.txt`. Canonical fetch: `antigenomics/vdjdb-db` GitHub releases. ⚠ **release matters** — the 2024-06 checkout has **zero** records for the AS/B27 motif (below); the 2025-12-29 dump has them (Yang 2022 landed in between). Do not answer "is X in VDJdb" from a stale dump |
| **Ankylosing spondylitis / SpA (`airr_ankspond`)** — AS/B27 motif campaign | HF dataset [`isalgo/airr_ankspond`](https://huggingface.co/datasets/isalgo/airr_ankspond) (cc-by-nc-nd-4.0, ~248 MB, TRB only). Two studies: **`old/`** = Komech et al. 2018, **`new/`** = Komech et al. 2022 | **`new/`** 68 `*.tsv.gz` (`count junction_aa v_call j_call locus`; sniffs `airr`) + `metadata.tsv` (`donor_id disease_status(as/hd) b27(pos/neg) proj sample_type time_point sample_name fraction`). **`old/`** 40 `*.txt.gz` in the **MiTCR/tcR dotted** dialect (`Read.count CDR3.nucleotide.sequence … V.end J.start D5.end D3.end VD.insertions`) + `metadata.tsv` (`id state(as/h) b27`) | `huggingface_hub.snapshot_download("isalgo/airr_ankspond", repo_type="dataset")`; read with `io.read` (`new/`→`airr`, `old/`→**`mitcr`**, the reader added for this cohort). Driver: `~/vcs/projects/2026-vdjtools-benchmark/bench/bm_ankspond.py` | **experimental** TCRβ repertoires. Komech et al., *Rheumatology* 2018;57(6):1097–1104, [10.1093/rheumatology/kex517](https://doi.org/10.1093/rheumatology/kex517) → `old/`; Komech et al., *Front Immunol* 2022;13:973243, [10.3389/fimmu.2022.973243](https://doi.org/10.3389/fimmu.2022.973243) → `new/`. **Gotchas**: `new/` is **amino-acid only** (`junction_nt` 100% null) — nt Pgen needs `old/`; `locus` is the literal string `"beta"` (vdjtools derives locus from `v_call` and ignores it); `"nan"` is the literal missing marker; `Koc.tsv.gz` has `Float64` counts where every other file is `Int64` (`io.read` normalizes; a naive `read_csv`+`concat` raises); the HF dataset **viewer is broken** for this repo (`load_dataset()` fails — use `snapshot_download`). **Design**: B27 is 26/27 confounded with AS (per-donor AS/B27+ 26, AS/B27− **1**, HD/B27+ 12, HD/B27− 21), so only the **B27-matched** contrast separates disease from carriage; `proj` marks 5 of the 12 HD/B27+ as foreign batches. **38 of 40 `old/` donors reappear in `new/`** (34 name-mapped agree on both phenotypes) — the two folders are **not** independent cohorts and must not be pooled |
| AS/B27 motif oracle (Yang 2022) | VDJdb **2025-12-29** dump, `/Users/mikesh/vcs/code/vdjdb-iedb-concordance/vdjdb_dump_2026/vdjdb.slim.txt` | TSV (as above) | filter `cdr3` ∈ the Komech 9; 7 records, all `TRBV9*01`/`TRBJ2-3*01`/`HLA-B*27:05`, score 3 | **curated** — Yang X, Garner LI, Zvyagin IV, …, Komech EA, …, Chudakov DM, …, Garcia KC. *Autoimmunity-associated T cell receptors recognize HLA-B\*27-bound peptides.* **Nature** 2022;612(7941):771–777, [10.1038/s41586-022-05501-7](https://doi.org/10.1038/s41586-022-05501-7), PMID 36477533 (verified via PubMed). Epitopes: `TRLALIAPK` (PRPF3, self), `GQVMVVAPR` (RNASEH2B, self), `LRVMMLAPF` (*E. coli* YEIH). Also supplies a **10th** family member absent from Komech 2018: `CASSVGTYSTDTQYF`. ⚠ **partly circular** — same group as Komech 2018, TCRs isolated from AS patients, plausibly these donors: it confirms restriction+epitope, it does **not** independently validate discovery |
| FMBA covid19 (TRA+TRB) — Phase 6b benchmark | aldan3 `/projects/fmba_covid` (clonotype tables `COV_V_usage_adjustment_v3/FMBA_functional/*.clonotypes.TRB.txt` + `data/*.clonotypes.{TRA,TRB}.pool.aa.table.txt`) ↔ private HF `isalgo/airr_covid19`; phenotype = `metadata_fmba_full.txt`/`desc_fmba_not_nan_hla.csv` (`COVID_status`, 4-digit `HLA-*`), join by 12-digit `id` | legacy VDJtools tables; TAB metadata | `~/vcs/projects/2026-vdjtools-benchmark/scripts/biomarker_bench.sbatch -- covid19` on aldan3 (reads `/projects/fmba_covid` directly) | **experimental** deep TCRα/β covid repertoires (Vlasova et al., *Genome Medicine* 2026;18:20). COVID-association + α-β co-occurrence benchmark |
| covid19 biomarker oracle (**canonical**) | VDJdb SARS-CoV-2 — same slim dump as the CMV row; on aldan3 `/projects/immunestatus/vdjdb/vdjdb-2025-07-30/vdjdb.slim.txt` | TSV (as above) | filter `species==HomoSapiens & antigen.species~"SARS-CoV-2"` (**both chains**: 3,796 TRA + 5,333 TRB records → 8,842 unique CDR3s); `bench_biomarker._vdjdb_antigen(path, "SARS-CoV-2")` | **curated** antigen-specific oracle (author decision 2026-07-16: prefer VDJdb over the study's published list — assay-grounded, chain-agnostic, symmetric with the CMV validation) |
| covid19 biomarker oracle (alternative, not used) | HF `isalgo/airr_covid19` → `covid_associated_clonotypes.csv` (`cdr3,cluster,has_covid_association,chain,v,j`) | CSV | `--oracle` flag (opt-in) | **curated** — the study's published COVID-associated clonotype/cluster list (Vlasova 2026). Superseded as the default validation target by VDJdb SARS-CoV-2 (above) |
| FMBA covid19_vacc (TRA+TRB) — Phase 6b benchmark | aldan3 `/projects/fmba_covid/vaccine/corr_func/*.clonotypes.{TRA,TRB}.txt` ↔ private HF `isalgo/airr_covid19_vacc`; phenotype = `vaccine/processed_metadata.tsv` / `vaccine/metadata.csv` (`timepoint` ∈ before_vaccination/20d_after_vaccination, `vaccine` ∈ GamCOVIDVac/CoviVac/EpiVacCorona) |  legacy VDJtools tables; CSV/TSV metadata | `~/vcs/projects/2026-vdjtools-benchmark/scripts/biomarker_bench.sbatch -- covid19_vacc` on aldan3 | **experimental** pre/post-vaccination TCRα/β repertoires (Vlasova 2026). Timepoint/vaccine association benchmark |

## Phase 14 — repertoire dynamics (paired tracking + neighbourhood enrichment)

Four longitudinal vaccination cohorts. Three of them had **zero references in this repo** before this
phase; none is on the cluster (`covid19_vacc` above is). All are per-sample gzipped TSV + a metadata
table — **the HF dataset viewer is broken or misleading for every one of them**, so use
`hf_hub_download`/`snapshot_download`, never `load_dataset()`.

| Dataset | Origin / fetch | Format | Structure | Provenance |
|---|---|---|---|---|
| **Yellow fever (`airr_yfv19`)** — the positive control | HF [`isalgo/airr_yfv19`](https://huggingface.co/datasets/isalgo/airr_yfv19); `huggingface_hub.snapshot_download("isalgo/airr_yfv19", repo_type="dataset")` | 42 × `<sample>.airr.tsv.gz` + **`metadata.txt`** (not `.tsv`; cols `file_name donor day replica`). Schema uses **`v_gene`/`d_gene`/`j_gene`** (not `*_call`) and has no `duplicate_frequency` — `io.read_airr` handles this as of Phase 14; before that V/J came back 100% null | 6 donors `P1 P2 Q1 Q2 S1 S2` = 3 monozygotic-twin pairs (inferable from the label prefix only — **no twin column**); days −1/0/7/15/45; **TRB only**; UMI counts (MIGEC+MiXCR per the paper; the repo documents nothing) | **experimental**. Pogorelyy MV et al., *PNAS* 2018;115(50):12704–12709, [10.1073/pnas.1809642115](https://doi.org/10.1073/pnas.1809642115), PMID 30459272. SRA PRJNA493983 |
| **TBEV (`airr_tbev_vac`)** — the only cohort shipping labels | HF [`isalgo/airr_tbev_vac`](https://huggingface.co/datasets/isalgo/airr_tbev_vac) | 287 × `samples/<id>.airr.tsv.gz` + `metadata.tsv` (34 cols; `v_call`/`j_call`, **read** counts not UMI) | **10 donor keys, not the README's 11** (`Vkh` aggregates the `VK_*`/`VKh_*` batches); ragged timepoint grid (`day` is the only unifying axis: `p0`/`0`/`0-1`/`0-2`/`1st` all → day 0); F1/F2 + a separate `y` replicate series; `_r2` = a second library; **12 `cell_subset`s** — a pair must match on subset, else it is a sorting artefact, not a timepoint effect; **class I + II HLA** for all 10 | **experimental**. Sycheva AL et al., *Front Immunol* 2022;13:970285, [10.3389/fimmu.2022.970285](https://doi.org/10.3389/fimmu.2022.970285), PMID 36091004. SRA PRJNA847436 |
| **Influenza (`airr_flu_vac`)** — the weak/negative control | HF [`isalgo/airr_flu_vac`](https://huggingface.co/datasets/isalgo/airr_flu_vac) | 17 × `samples/<id>.airr.tsv.gz` + `metadata.tsv` (19 cols, **CRLF**; `v_call`/`j_call` **gene-level only**, no allele) | **ONE donor (`YB`)** — no cohort-level analysis is possible; seasons 13/14/16 are the longitudinal axis; timepoints `pre0`/`0`/`5`/`12`/`45` (day −14/0/5/12/45); F1/F2 replicates in seasons 14+16 only; **UMI-corrected** `duplicate_count`; class I HLA (one genotype) | **experimental**. Sycheva AL et al., *Vaccine* 2018;36(12):1599–1605, [10.1016/j.vaccine.2018.02.027](https://doi.org/10.1016/j.vaccine.2018.02.027), PMID 29454515. SRA SRP111073. ⚠ **paywalled, no PMC copy** — the abstract reports "several" new memory clonotypes peaking day 45 and gives no method or count, so **there is no verified numeric ground truth** for this cohort |
| **YF-specific TCR reference** (the enrichment validation target) | VDJdb via `vdjmatch.db.load(asset="slim", gene="TRB")`, filtered `antigen.epitope == "LLWNGPMAV"`. NB `vdjmatch.db.fetch_hf` **defaults to `repo="isalgo/airr_benchmark"`** — an upstream provenance dependency this file did not previously record | TSV (`gene cdr3 antigen.epitope v.segm j.segm mhc.a …`) | **409 unique TRB CDR3aa**, A\*02:01/NS4b-restricted | **curated, and genuinely held-out**: every record is `tetramer-sort` or `dextramer-sort` (VDJdb `method.identification`; **zero** expansion- or frequency-derived entries), from five independent studies (Lee 2017 PMID 28103239; Afik 2017 PMID 28934479; Bovay 2018 PMID 28975614; Rius 2018 PMID 29483360; vdjdb-db#193) — **no Pogorelyy paper, no twin cohort**. ⚠ `mhc.a` resolution is inconsistent (486 `HLA-A*02:01`, 188 `HLA-A*02`) — normalise before filtering; `vdjdb.score` is 0 for 391/409, so `min_score>=1` leaves 18 and kills the test |

⚠ **`airr_benchmark/alice/yf/{P1..S2}_{d0,d15}.tsv.gz` is NOT a reference set.** It is the ALICE paper's
benchmark *repertoires* — the same 6 donors and the same two timepoints as `airr_yfv19`, i.e. the thing
under test. Using it as labels would be circular. The reference is the VDJdb row above.

⚠ **`airr_yfv19` traps, all verified against the actual metadata rather than the paper:**
`day = -1` is a **mislabel** for the paper's day −7 (a 7× error on that interval for anything
time-weighted). **"Two replicates per timepoint" is false** — only 12 of 30 (donor, day) cells have both
(F1 ×30, F2 ×12; only S2 has the full design). But **every donor has both day-0 replicates**, which is
exactly what the day0→day15 analysis needs. **No ground truth ships**: no expanded/responding lists, no
dextramer/IFNγ/activated fractions, no HLA (the paper's 600–1700 responding TCRβ/donor and the 395/773
donor-S1 validation live in its SI + github.com/mptouzel/pogorelyy_et_al_2018). ~4% of rows are
**non-functional** (`*`/`_` in `junction_aa`) and are retained. The viewer "works" but concatenates all
42 files into one flat 31.5M-row table **with no donor/day/replica column** — it silently destroys the
longitudinal design.

⚠ **`airr_tbev_vac`**: gate on `reads`, **not** on `library=main` — `IZ_reactive_IFNgpos` is `main` at
2,384 reads. The antigen-reactive sorted fractions (IFNγ±, IL2±) are the only shipped labels anywhere in
the four cohorts and are also the **shallowest** samples (min 439 reads).

### Phase 6b benchmark results (Aldan-3, v2.7.0, 16 cores) — **computed**, not experimental

**The scripts moved out on 2026-07-17** — they live in `~/vcs/projects/2026-vdjtools-benchmark` (`bench/` + `scripts/`); the numbers of record stay here.

All three via `~/vcs/projects/2026-vdjtools-benchmark/scripts/biomarker_bench.sbatch`; `key=(junction_aa,v_call,j_call)`, BH-FDR q<0.05,
`alternative="greater"`. Concordance = how many of the top-100 Fisher hits each test also flags.

**Arms are the ANALYSED counts** (`association()` inner-joins the design to the cohort — the design
size is *not* the tested size; see `run_association_suite`, which now reports and warns on this).

| Cohort | Analysed arms | min_inc | Features tested | Significant | Validation | χ² | perm | BF | Time / peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| covid19 (COVID vs healthy) | **502+ / 34−** ⚠ 15:1 | 10 | 52,528 | **0** | none to validate (OR=nan) — 34 controls cannot power a genome-wide screen | — | — | — | 142 s / 6.4 GB |
| **hip (Emerson CMV)** — *the association baseline* | 340+ / 421− | 8 | **1,366,592** | **70** | **VDJdb-CMV OR=24.5, p=5.5e-11** (10/70 sig CDR3s known; 4,551/671,879 tested are members); CMH (CMV\|HLA-A\*02) → 40 significant | 39/100 | 0/100 | **100/100** | 2,536 s / 99.7 GB |
| ~~covid19_vacc (timepoint)~~ **WITHDRAWN** | 541+ / 541− | 5 | 1,390,129 | ~~269~~ | **not quotable — see below** | — | — | — | 3,678 s / 43.7 GB |

⚠ **covid19 association is an honest negative: 0 significant on 502 COVID+ vs 34 healthy.** The FMBA
cohort is overwhelmingly COVID+ and only 34 healthy subjects have both metadata and a repertoire — too
few to power a genome-wide screen. **Earlier revisions of this file reported 44,125 significant with a
VDJdb-SARS-CoV-2 enrichment of OR=2.35 (p=3.5e-13). Both are withdrawn**: they were produced by the
phantom-negative bug (`association()` derived n_pos/n_neg from the *design* frame, so the 640 labelled
-but-unsequenced subjects in the 1,212-row FMBA sheet voted "feature absent" for every feature — 438 of
472 negatives were phantoms). Fixed in v2.7.0; on fixed code the count is **0**. hip (340/421) and
covid19_vacc (541/541) never had designs that outran their cohorts and are unaffected — **hip is the
association baseline**, covid19 is not.

⚠ **covid19_vacc (timepoint) is withdrawn for a different and unfixed reason: the design.** The 269
significant features were produced by an *incidence* test — "is this clonotype carried by more
after-samples than before-samples?" — asked of a **paired** design. Two things follow, and neither is a
bug in the code that produced them:

1. **A 29% richness gap survives depth equalisation.** At identical read counts the post-vaccination
   repertoires carry ~29% more unique clonotypes (pre-downsample 12,854 vs 9,936; post-downsample
   11,990 vs 9,295 — *unchanged*). A richer repertoire contains any given clonotype more often, so an
   incidence test reports that richness, not the vaccine. Per-pair read-downsampling does not touch it.
2. **The gap is compositional, not just sampling.** Because frequencies sum to 1, ~2,900 extra clones at
   ~8–10 reads each take ε ≈ 3–6% of the mass, so *every* shared clonotype is deflated by (1−ε) in the
   richer sample. H₀ is false for all of them, worst where `N_eff·f·ε²` is large — i.e. at the top of the
   distribution, reproducibly, in all 541 donors.

The paired design supports a **within-donor frequency** question instead (`vdjtools.dynamics`), which
conditions the unknown frequency away exactly (thesis Eq. 2.4). That addresses (1). It does **not**
address (2), which is orthogonal to the sampling noise `N_eff` models — so the arm is re-run to *measure
and report* the compositional bias, not to produce a corrected count. Note edgeR's TMM factor is an
estimate of ε, and edgeR is already a committed comparator.

**Donor key**: the donor is `id`; `name` is the 12-digit **per-sample** id (donor × timepoint) and
`full_id` = `<id>_<pre|post>_<vaccine>`. Using `name` as the donor makes every sample its own donor.
Earlier revisions of this row said "join by 12-digit `id`", conflating the two.

Depth does **not** confound the association: repertoire size is unrelated to COVID status (medians
14,146 vs 13,211, ratio 1.07×, Mann-Whitney **p=0.64**, `~/vcs/projects/2026-vdjtools-benchmark/bench/assoc_depth.py`), so no depth
conditioning is applied there — unlike co-occurrence, where depth is the whole problem.

α-β **co-occurrence** (covid19, `min_cooccurrence=3`, `min_incidence_frac=0.03`, `evalue=True`,
**`depth_strata=10` default**): 481,538 candidate pairs → **502 significant** (q<0.05); θ median
**3.52**, max **9.19**; median depth-conditioned `or_mh` **8.97**; 17 s. Validated against **VDJdb α-β
complexes (any antigen**, 20,169 pairs — pairing is a receptor property, so the oracle is deliberately
*not* antigen-restricted): **3/502 significant pairs are known VDJdb pairs vs 6/446,224 tested →
enrichment OR=893.2, p=2.8e-08**.

**The depth correction improved precision without costing recall**: it removed 82% of the pooled hits
(2,802 → 502) while retaining **all 3** VDJdb-validated true pairs, so the enrichment rose 5.6×
(OR 158.4 → 893.2; p 4.9e-06 → 2.8e-08). Removing artifact, not signal. Honest caveat: only 6 known
pairs were testable, so this rests on a small oracle overlap (Fisher is exact, so the p is valid).
⚠ candidate features capped at `max_features=2000` (warned at runtime), i.e. top-incidence pairs only.

#### Co-occurrence confounding — measured, and corrected by default (`appendix/{depth_gate,theta_ceiling,accept_gate}.py`)

Cross-subject co-occurrence is confounded by per-subject **repertoire depth** and **shared HLA**
(see the caveat in `docs/usage.rst`). Depth was quantified, then **fixed**; HLA remains a caveat.

**Depth.** Subjects span **1 → 90,174** unique clonotypes (CV=0.899), inducing a lift
`θ_depth = 1+CV² = 1.809` for rare clonotypes with no biology whatsoever. `~/vcs/projects/2026-vdjtools-benchmark/bench/cooccurrence_fpr.py`
(seed 20260716; 2,000 null pairs/config; depth lognormal at the cohort's CV) measures what that does
to a **pooled** test — and `max_features=2000` keeps the top features *by incidence*, so the tool
operates in the worst column:

| false-positive rate @ nominal p<0.05 (calibrated = 0.05) | 2% incidence | 3% | **11%** |
|---|---|---|---|
| pooled Fisher (the pre-2.7.0 default) | 0.052 | 0.073 | **0.456** |
| **CMH over depth strata (the 2.7.0 default)** | **0.023** | **0.022** | **0.036** |

The **depth-weighted / "x of X rearrangements"** variant was evaluated and **rejected as not a
well-defined test**, not on its rate: a hypergeometric is a statement about drawing discrete units,
so feeding it sums of depths silently asserts ΣS≈7e6 exchangeable pseudo-observations in place of 552
subjects. Its "FPR" is then an artifact of the variant written — the formulation in
`cooccurrence_fpr.py` degenerates conservative (→0.000), while one rescaling the weights to keep *n*
fixed degenerates anticonservative (→0.9). Unusable in either direction. (An earlier revision of this
file quoted 0.302/0.461/0.887 as measured; that number is not reproducible and is withdrawn.)

**⚠ Correction.** An earlier revision of this file argued the signal was *"not a depth artifact —
94.1% of significant pairs exceed θ_depth=1.809"*. **That reasoning is invalid**: `θ_depth=1+CV²`
is the null **mean**, not a tail quantile, so ~half of *pure-null* pairs exceed it by construction
and conditioning on significance selects the high-θ ones. ~94% is what a null screen *should*
produce. The claim is withdrawn; the table above measures the thing that actually matters.

**Effect of the fix on the real cohort** (`~/vcs/projects/2026-vdjtools-benchmark/bench/accept_gate.py`, all 552 subjects; **one
deterministic run** — both arms share the identical candidate set, so this is a controlled A/B in
which only the conditioning differs):

| covid19 α-β co-occurrence | tested | significant (q<0.05) | θ median / max | median `or_mh` | time |
|---|---|---|---|---|---|
| `depth_strata=0` (pooled, uncorrected) | 481,579 | **2,804** | 3.07 / 9.19 | — | 16 s |
| **`depth_strata=10` (default)** | 481,579 | **641** | 3.38 / 9.19 | **8.30** | 10 s |
| ≥1000-clonotype floor + pooled (diagnostic) | 733,178 | 754 | 3.47 / 7.57 | — | — |

Depth conditioning removes **77%** of the pooled hits *without discarding a subject*, and the
survivors have a higher θ and a strong depth-conditioned `or_mh`. **641 is the defensible figure.**
The depth floor is now a diagnostic, not required preprocessing.

⚠ **Earlier revisions of this file reported 2,805 / 502 / "82%" and three different tested counts
(481,538 / 481,501 / 481,587) for one configuration.** Those predate two fixes and are withdrawn:
candidate selection was nondeterministic (an incidence tie-break, so `max_features` cut a different
set every call — the arms were never comparable), and the CMH branch ignored `alternative`, making
the "corrected" arm a *two-sided* test whose hits included anti-correlated pairs. The table above is
a single post-fix run with a stable tested count in both arms.

**HLA is the remaining binding constraint — untouched by this.** Shared restriction by an allele of
carrier frequency *f* cannot induce a lift above `(1+CV²)/f`; only the extreme tail (θ=9.19) clears
every allele carried by >19.7% of the 466 typed subjects (A\*02:01, A\*03:01, A\*01:01, A\*24:02,
B\*07:02). Ancestry, batch and shared exposure are likewise unfixed. **None of this licenses a
physical-pairing claim** — a cross-subject α-β pair is pairSEQ's *definition of a false positive*
(Howie 2015). TODO v2.8: a degree-preserving (fixed-fixed) permutation null — the truest null
measured (lift 0.98–1.01) but ~120 core-days at 1e-6 resolution, hence not the production default.

**Three scale effects** (reproduced independently across cohorts — properties of genome-wide incidence
testing, not defects): (1) BH-FDR is severe at ~10⁶ features — Emerson-2017 itself thresholds at a
*nominal* P<1e-4 (FDR≈0.14), so 70 CMV hits at q<0.05 is consistent with the paper; (2) **permutation has
a 1/`n_perm` p-floor** — at `n_perm=1000` the smallest p is 1e-3, and BH over 1.37M features needs
≈1e-3×1.37M/70≈19 ≫ 0.05, hence 0/100 on both ~1.4M-feature cohorts but 100/100 on covid19's 52k set;
(3) Yates-corrected **χ² is asymptotic** and conservative on the sparse tables at low `min_incidence`
(39/100 → 85/100 → 100/100 as counts grow). **Fisher and the Beta-Binomial BF are the reliable arbiters
at genome-wide scale; χ² and permutation are scale-dependent.**

## Phase 5 — preprocess (VJ batch-correction validation)

| Dataset | Origin | Format | How to obtain | Provenance |
|---|---|---|---|---|
| FMBA covid Cohort I (TCRβ) | aldan3 HPC project `fmba_covid`: clonotype tables `COV_V_usage_adjustment_v3/FMBA_functional/*.clonotypes.TRB.txt`; raw usage `COV_V_usage_adjustment_v3/V_usage_FMBA.tsv`; paper-adjusted `…_adjusted.tsv`; batch labels = `unprocessed_fmba_metadata.csv` col `folder` (9 NovaSeq runs = the paper's 9 batches), join by 12-digit id in col `name` | legacy VDJtools tables (`count freq cdr3nt cdr3aa v d j VEnd…`); usage matrices = sample × V-gene/family TSV | `aldan3 pull /projects/fmba_covid/COV_V_usage_adjustment_v3/…` (needs the `aldan3` client + `fmba_covid` group). Validation driver: `~/vcs/projects/2026-vdjtools-benchmark/bench/validate_batch_covid.py` | **experimental** — deep TCRβ repertoires (median ~3.3M reads/sample) from Vlasova et al., *Genome Medicine* 2026;18:20 (10.1186/s13073-025-01589-4). Used to validate `preprocess.correct_vj_usage(transform="sigmoid")` + `apply_vj_correction`: batch η² 0.109→0.002, grand mean preserved, reads preserved |

## Phase 10 — BCR SHM & lineage (test data, planned; not yet fetched)

Verified via PubMed (2026-07); all four are owner (Shugay)-co-authored BCR datasets. Data paths are
**TBD** — fill the `How to obtain` details when each dataset is actually sourced; do not fabricate a
path before then. All are **experimental** repertoire data.

| Dataset | Reference (verified) | DOI | Relevance |
|---|---|---|---|
| Longitudinal memory-B / ASC BCR repertoires | Mikelov et al., *eLife* 2022;11:e79254 | [10.7554/eLife.79254](https://doi.org/10.7554/eLife.79254) | UMI-tagged BCR-seq + clonal-lineage/phylogeny — primary Phase-9 (MiGEC) + Phase-10 (lineage) target |
| TCGA tumor repertoires from RNA-seq | Bolotin et al., *Nat. Biotechnol.* 2017;35(10):908–911 | [10.1038/nbt.3979](https://doi.org/10.1038/nbt.3979) | RNA-seq-derived TCR/BCR (MiXCR method); TCGA IGH |
| CD27-dull/bright memory-B VH repertoires | Grimsholm et al., *Cell Rep.* 2020;30(9):2963–2977.e6 | [10.1016/j.celrep.2020.02.022](https://doi.org/10.1016/j.celrep.2020.02.022) | memory-B subsets, VH usage + SHM frequency (owner wrote "CD20-dull"; the paper is **CD27**-dull) |
| CVID peripheral B-cell selection | Grimsholm et al., *Cell Rep.* 2023;42(5):112446 | [10.1016/j.celrep.2023.112446](https://doi.org/10.1016/j.celrep.2023.112446) | CVID Ig-seq, peripheral B-cell selection |

**Flag (unresolved):** the owner's "Mikelov *allergy* paper" literally matches Mikelov et al.,
*Nat. Immunol.* 2025;26(12):2328–2342 ([10.1038/s41590-025-02323-3](https://doi.org/10.1038/s41590-025-02323-3),
peanut oral-immunotherapy) — but that study is **TCR / single-cell** and does **not** list Shugay
as an author, so the Shugay-co-authored **eLife 2022 BCR** paper is recorded above instead. Confirm
which Mikelov dataset is intended before use. (References verified via PubMed; DOIs copied verbatim.)

## Bundled precomputed models (shipped in the wheel)

Ship under `python/vdjtools/model/_bundled/<source>/<LOCUS>/` (parquet marginals + `manifest.json`);
loaded with `vdjtools.model.load_bundled(locus, source)`. ~0.4 MB total (30–150 KB/model).

| Model set | Origin | How to rebuild | Provenance |
|---|---|---|---|
| `olga` (7 loci) | `from_olga` on the OLGA default models above | `python appendix/build_bundled_models.py` (the OLGA part is inline in the precompute) | **derived** — OLGA generative parameters converted to the polars schema (keeps OLGA germline; exact-Pgen bootstrap, single-D) |
| `learned` (7 loci) | native EM (`infer_native`) on real out-of-frame HF reads (2 000 unique clonotypes/locus, 12 fixed iters, arda V/J[/D]-masked) seeded from the `olga` model | `python appendix/build_bundled_models.py` | **computed from experimental data** — real-repertoire gene-usage/trim/insertion marginals. D-bearing loci (IGH/TRD/TRB) carry an **arda-anchored tandem-D** event: a read may be `n_D=2` only where arda called a second D (`d2_call`), which counters the tandem-vs-long-insertion identifiability that inflates unregularized D-D EM (TRB 0.28→**0.00**, TRD 0.18→**0.006**, IGH **0.009**). `EM_SINGLE_D=1` / `ND_PRIOR` are alternative regularizers. |

## Golden fixtures (tests)

| Dataset | Origin | Use |
|---|---|---|
| IGoR model files (human/mouse, all loci) | `IGoR-models/` (local) | model-loader + Pgen oracle fixtures; canonical `model_parms`/`model_marginals`/anchors format |
| OLGA `default_models/` (7 human loci) | **vendored in-repo** at `tests/python/fixtures/olga/default_models/` — copied from `antigenomics/mirpy` (`main`/`legacy-v2`, `mir/resources/olga/default_models/`; byte-identical to OLGA's published models, Sethna et al.). Re-fetch: `git -C ../mirpy archive main mir/resources/olga/default_models \| tar -x`. Mouse loci not vendored (unused). | model-loader + native-vs-OLGA Pgen oracle (`from_olga`). Tests resolve them here (override with `VDJTOOLS_OLGA_MODELS`); the `olga` pip package (`[oracle]` / test extra) supplies the *runtime* Pgen comparison. Makes the suite self-contained — no external mirpy checkout needed. |
| Legacy input-format samples | `tests/python/fixtures/legacy/*.txt.gz` — copied from legacy vdjtools `src/test/resources/samples/` (on the `legacy-1.x` branch; re-fetch with `git show legacy-1.x:src/test/resources/samples/<f>.txt.gz`) | format-conversion conformance for `vdjtools.io.convert` (MiXcr v1/2+v3, MiGec, ImmunoSeq v1/v2, ImgtHighVQuest, Vidjil, RTCR). Oracles (input row → canonical values) derived from the legacy Groovy parsers; see `tests/python/test_convert.py` |
| MiTCR/tcR sample | `tests/python/fixtures/legacy/mitcr.txt.gz` — **real data slice** (not synthetic): header + first 8 rows copied verbatim from `isalgo/airr_ankspond` `old/Luk4y.txt.gz` (the smallest `old/` file), the MiTCR output of Komech 2018. Re-cut: `gzip.open(<snapshot>/old/Luk4y.txt.gz).readlines()[:9]` | conformance for `read_mitcr` + `sniff_format` → `"mitcr"` (`tests/python/test_convert.py`). Oracle rows verified against the raw file, incl. the ambiguous `D.gene` = `"TRBD1, TRBD2"` → `extract_vdj` keeps `TRBD1` |
| TRUST4 report sample | `tests/python/fixtures/legacy/trust4.txt.gz` — **synthetic** (not from legacy vdjtools), hand-built from the literal TRUST4 `*_report.tsv` header (`#count frequency CDR3nt CDR3aa V D J C cid cid_full_length`) in [`liulab-dfci/TRUST4`](https://github.com/liulab-dfci/TRUST4) `trust-simplerep.pl`; CDR3 nt/aa reuse verified junctions from the MiXcr/RTCR fixtures | conformance for `read_trust4` (`tests/python/test_convert.py`), incl. `partial`/gene-less row skipping. arda AIRR output is exercised inline (arda writes standard AIRR; `read_arda` delegates to `read_airr`) |

## Carried-over resource data (from legacy vdjtools, `legacy-1.x` branch)

| File | Use |
|---|---|
| `resources/profile/aa_property_table.txt` | amino-acid physicochemical properties (Phase 3) |
| `resources/profile/cdr3contact.txt` | CDR3 contact-probability estimate (Phase 3) |
