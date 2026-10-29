"""Build the EM-learned bundled models from real HuggingFace non-functional reads, all 7 loci.  2026-07-12

For each locus: fetch the FULL non-functional read set (HF isalgo/airr_model_read), arda-map to
unique clonotypes = (v_call, j_call, junction), and run native EM (D-D by default on the D-bearing
loci IGH/TRD/TRB) seeded from the ARDA germline scaffold. Saves to
python/vdjtools/model/_bundled/learned/<locus>/ (parquet marginals + manifest.json + training.json).

The pipeline itself now lives in the library — `vdjtools.model.data.build_model` / `build_all`,
also reachable as `vdjtools model build` — so this script only supplies the two things specific to
*rebuilding what ships in the wheel*: the organism/locus selection and the destination directory. Anything
that merely wants to train a model should call the library directly.

Non-functional = out-of-frame OR stop-codon, and BOTH are used. The only property the generative
model needs is that the rearrangement escaped selection, and both halves have. Keeping only the
out-of-frame half would condition the training set on junction length mod 3 — a bias the
insertion-length model would then happily learn.

No cap, no subsampling: every clonotype surviving the germline filter goes into EM.

Convergence-based EM: stop when the relative log-likelihood improvement falls below EM_TOL
(whole-model signal — the old V-usage criterion settled in ~2 iters while trims/insertions/n_d were
still moving, so it was disabled; log-lik is monotone now that the fixes keep the scoreable-read set
stable). EM_ITERS is a generous safety cap, not the target.

Reproduce (needs arda/mmseqs2 + HF access):  python appendix/build_bundled_models.py

For arda-native models set BUNDLE_SOURCE=arda and BUILD_GROUPS=human,mouse.
LOCI selects affected loci. RAW_FASTQ_DIR points to a local raw-input mirror with
<group>/<locus>.nonfunctional.fq.gz; those reads are freshly annotated, without a cap.
"""
import os
import time
from pathlib import Path

import polars as pl

from vdjtools.model import data

SOURCE = os.environ.get("BUNDLE_SOURCE", "learned")
DEST = Path("python/vdjtools/model/_bundled") / SOURCE
# NB the previous default pointed into ANOTHER SESSION's scratchpad and the build silently
# reused a truncated arda cache from it -- that is how the shipped TRB model came to be
# trained on 870 clonotypes when 32,562 were available. arda.annotate_reads does its own
# caching under out_dir; there is no second cache layer here any more.
WORK = Path(os.environ.get("EM_WORK", "/tmp/em_work"))
LOCI = ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL")
MOUSE_LOCI = ("TRA", "TRB", "IGH", "IGK", "IGL")
ITERS = int(os.environ.get("EM_ITERS", "15"))          # generous safety cap
EM_TOL = float(os.environ.get("EM_TOL", "1e-4"))       # stop at relative log-lik improvement < this
WORKERS = int(os.environ.get("EM_WORKERS", "0")) or None
THREADS = int(os.environ.get("EM_THREADS", "0"))


def training_clones(locus: str, group: str = "human"):
    """The clonotypes EM trains on, from the pre-annotated cache when one is staged.

    EM_DATA_DIR: read pre-annotated slim parquet (v/j/d/junction/...) instead of fetching from
    HuggingFace + re-annotating with arda. This is the Aldan-3 path -- the compute nodes have no
    mmseqs2 and no outbound HTTPS, so annotation happens once on a workstation and the columns EM
    needs are staged as compact parquet (~20 MB for all 7 loci vs 1.6 GB of raw arda TSV).
    """
    raw_dir = os.environ.get("RAW_FASTQ_DIR")
    data_dir = os.environ.get("EM_DATA_DIR")
    if raw_dir and data_dir:
        raise ValueError("select raw FASTQ inputs or pre-annotated inputs, not both")
    if raw_dir:
        reads = data.annotate_reads(
            Path(raw_dir) / group / f"{locus}.nonfunctional.fq.gz",
            out_dir=WORK / group / locus, prefix=f"{group}_{locus}",
            organism=data.ORGANISM[group], threads=THREADS, verbose=True)
        return data.unique_clonotypes(reads, gene_level=True)
    if not data_dir:
        return None
    return data.unique_clonotypes(
        pl.read_parquet(Path(data_dir) / f"{group}_{locus}_nonfunctional.parquet"), gene_level=True)


def main():
    if SOURCE not in ("learned", "arda"):
        raise ValueError("BUNDLE_SOURCE must be learned or arda")
    groups = os.environ.get("BUILD_GROUPS", "human,mouse").split(",")
    if any(group not in ("human", "mouse") for group in groups):
        raise ValueError("BUILD_GROUPS must select human and/or mouse")
    want = os.environ.get("LOCI")
    loci = tuple(want.split(",")) if want else LOCI
    if set(loci)-set(LOCI):
        raise ValueError("LOCI must select supported receptor loci")
    if want and "mouse" in groups and set(loci)-set(MOUSE_LOCI):
        raise ValueError("mouse bundled training targets are TRA/TRB/IGH/IGK/IGL")
    DEST.mkdir(parents=True, exist_ok=True)
    # EM_SINGLE_D=1 forces strict single-D; ND_PRIOR adds a Dirichlet single-D pseudocount.
    # GENE_PRIOR keeps every germline-functional V/J allele reachable -- P(V)=0 is an EM absorbing
    # state, so unregularized, human TRB kept 30 of OLGA's 57 V genes having SEEN 54 in the data.
    single_d = os.environ.get("EM_SINGLE_D") == "1"
    nd_prior = float(os.environ.get("ND_PRIOR", "0"))
    gene_prior = float(os.environ.get("GENE_PRIOR", "1.0"))

    started = time.perf_counter()
    rows = []
    for group in groups:
        for locus in loci:
            if group == "mouse" and locus not in MOUSE_LOCI:
                continue
            key = f"{group}_{locus}" if SOURCE == "arda" or group != "human" else locus
            print(f"[{key}] building {SOURCE} model ...", flush=True)
            model, rep, stats = data.build_model(
                locus, group=group,
                template=None,  # data.build_model owns the ARDA germline scaffold
                clones=training_clones(locus, group), threads=THREADS, verbose=True,
                work_dir=WORK, iters=ITERS, tol=EM_TOL, single_d=single_d,
                nd_prior=nd_prior, gene_prior=gene_prior)
            model.save(DEST / key)
            rows.append(stats)
            print(f"  [{locus}] {stats['n_used']}/{stats['n_clonotypes']} clonotypes, "
                  f"{stats['iters']} iters, {stats['seconds']:.0f}s, P(n_D=2)={stats['p_nd2']}, "
                  f"TRAINING LL {stats['loglik_first']:.0f}->{stats['loglik_last']:.0f}", flush=True)

    if rows:
        print(f"\n{pl.DataFrame(rows)}")
    print(f"\ntotal {time.perf_counter() - started:.0f}s for {len(rows)} models")


if __name__ == "__main__":
    main()
