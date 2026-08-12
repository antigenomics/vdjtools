"""Build the EM-learned bundled models from real HuggingFace non-functional reads, all 7 loci.  2026-07-12

For each locus: fetch the FULL non-functional read set (HF isalgo/airr_model_read), arda-map to
unique clonotypes = (v_call, j_call, junction), and run native EM (D-D by default on the D-bearing
loci IGH/TRD/TRB) seeded from the OLGA bootstrap model. Saves to
python/vdjtools/model/_bundled/learned/<locus>/ (parquet marginals + manifest.json + training.json).

The pipeline itself now lives in the library — `vdjtools.model.data.build_model` / `build_all`,
also reachable as `vdjtools model build` — so this script only supplies the two things specific to
*rebuilding what ships in the wheel*: the OLGA seed models, and the destination directory. Anything
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
"""
import os
import time
from pathlib import Path

import polars as pl

from vdjtools.model import data, from_olga

# The OLGA models shipped in THIS repo, not pip olga's: pip ships only 5 human loci (no
# TRG/TRD) plus mouse, while tests/python/fixtures/olga/default_models carries all 7 human loci.
# The TRG/TRD marginals originate from mirpy's legacy-v2 branch (commit aeccd75) and are verified
# byte-identical to what the bundled parquet were built from; olga-pip scores with them fine, so
# they are a real oracle for those two loci, which pip alone cannot be.
_REPO_OLGA = Path(__file__).resolve().parent.parent / "tests" / "python" / "fixtures" / "olga" / "default_models"
OLGA = Path(os.environ.get("VDJTOOLS_OLGA_MODELS", str(_REPO_OLGA)))
DEST = Path("python/vdjtools/model/_bundled/learned")
# NB the previous default pointed into ANOTHER SESSION's scratchpad and the build silently
# reused a truncated arda cache from it -- that is how the shipped TRB model came to be
# trained on 870 clonotypes when 32,562 were available. arda.annotate_reads does its own
# caching under out_dir; there is no second cache layer here any more.
WORK = Path(os.environ.get("EM_WORK", "/tmp/em_work"))
LOCI = {"TRA": "human_T_alpha", "TRB": "human_T_beta", "TRG": "human_T_gamma", "TRD": "human_T_delta",
        "IGH": "human_B_heavy", "IGK": "human_B_kappa", "IGL": "human_B_lambda"}
ITERS = int(os.environ.get("EM_ITERS", "15"))          # generous safety cap
EM_TOL = float(os.environ.get("EM_TOL", "1e-4"))       # stop at relative log-lik improvement < this
WORKERS = int(os.environ.get("EM_WORKERS", "0")) or None


def seed_model(locus: str, name: str):
    """The OLGA bootstrap model this locus's EM starts from.

    derive_orf=True: reconstruct the CDR3 germline for ORF alleles OLGA leaves empty (TRBV23-1 is
    8.6% of real TRB), so EM can learn their usage instead of pinning them to P(V)=0. Safe here --
    this is an arda-native learned model, not the exact-OLGA Pgen oracle (which keeps derive_orf off).
    """
    return from_olga(OLGA / name, locus=locus, derive_orf=True)


def training_clones(locus: str):
    """The clonotypes EM trains on, from the pre-annotated cache when one is staged.

    EM_DATA_DIR: read pre-annotated slim parquet (v/j/d/junction/...) instead of fetching from
    HuggingFace + re-annotating with arda. This is the Aldan-3 path -- the compute nodes have no
    mmseqs2 and no outbound HTTPS, so annotation happens once on a workstation and the columns EM
    needs are staged as compact parquet (~20 MB for all 7 loci vs 1.6 GB of raw arda TSV).
    """
    data_dir = os.environ.get("EM_DATA_DIR")
    if not data_dir:
        return None
    return data.unique_clonotypes(
        pl.read_parquet(Path(data_dir) / f"human_{locus}_nonfunctional.parquet"))


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    want = os.environ.get("LOCI")
    loci = {k: LOCI[k] for k in want.split(",")} if want else LOCI
    # EM_SINGLE_D=1 forces strict single-D; ND_PRIOR adds a Dirichlet single-D pseudocount.
    # GENE_PRIOR keeps every germline-functional V/J allele reachable -- P(V)=0 is an EM absorbing
    # state, so unregularized, human TRB kept 30 of OLGA's 57 V genes having SEEN 54 in the data.
    single_d = os.environ.get("EM_SINGLE_D") == "1"
    nd_prior = float(os.environ.get("ND_PRIOR", "0"))
    gene_prior = float(os.environ.get("GENE_PRIOR", "1.0"))

    started = time.perf_counter()
    rows = []
    for locus, name in loci.items():
        print(f"[{locus}] building learned model ...", flush=True)
        try:
            model, rep, stats = data.build_model(
                locus, template=seed_model(locus, name), clones=training_clones(locus),
                work_dir=WORK, iters=ITERS, tol=EM_TOL, single_d=single_d,
                nd_prior=nd_prior, gene_prior=gene_prior)
            model.save(DEST / locus)
            rows.append(stats)
            print(f"  [{locus}] {stats['n_used']}/{stats['n_clonotypes']} clonotypes, "
                  f"{stats['iters']} iters, {stats['seconds']:.0f}s, P(n_D=2)={stats['p_nd2']}, "
                  f"TRAINING LL {stats['loglik_first']:.0f}->{stats['loglik_last']:.0f}", flush=True)
        except Exception as e:  # noqa: BLE001 — report and continue
            print(f"  [{locus}] FAILED {type(e).__name__}: {e}", flush=True)

    if rows:
        print(f"\n{pl.DataFrame(rows)}")
    print(f"\ntotal {time.perf_counter() - started:.0f}s for {len(rows)}/{len(loci)} loci")


if __name__ == "__main__":
    main()
