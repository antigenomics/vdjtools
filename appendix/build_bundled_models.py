"""Build the EM-learned bundled models from real HuggingFace non-functional reads, all 7 loci.  2026-07-12

For each locus: fetch the FULL non-functional read set (HF isalgo/airr_model_read), arda-map to
unique clonotypes = (v_call, j_call, junction), and run native EM (D-D by default on the D-bearing
loci IGH/TRD/TRB) seeded from the OLGA bootstrap model. Saves to
python/vdjtools/model/_bundled/learned/<locus>/ (parquet marginals + manifest.json).

Non-functional = out-of-frame OR stop-codon; BOTH are used. The only property the generative model
needs is that the rearrangement escaped selection, and both halves have. Keeping only the
out-of-frame half would condition the training set on junction length mod 3 — a bias the
insertion-length model would then happily learn.

No cap, no subsampling: every clonotype surviving the germline filter goes into EM.

Fixed-iteration EM (tol=0): arda masks pin V, so the V-usage convergence check would stop after ~2
iters before the trims/insertions/n_d converge — running a fixed budget converges them properly.

Reproduce (needs [model] extra + arda/mmseqs2 + HF access):  python appendix/build_bundled_models.py
"""
import os
import time
from pathlib import Path

import polars as pl

from vdjtools.model import data, from_olga
from vdjtools.model.infer import _gene_to_alleles, gene_masks, infer_native

import olga as _olga
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
ITERS = int(os.environ.get("EM_ITERS", "12"))


def build(locus: str, name: str) -> dict:
    base = from_olga(OLGA / name, locus=locus)
    dset = set(base.genomic["genes_d"]["d_allele"].to_list()) if base.chain_type == "VDJ" else set()
    vgenes = set(_gene_to_alleles(base, "v"))
    jgenes = set(_gene_to_alleles(base, "j"))

    # EM_DATA_DIR: read pre-annotated slim parquet (v/j/d/junction/...) instead of fetching from
    # HuggingFace + re-annotating with arda. This is the Aldan-3 path -- the compute nodes have no
    # mmseqs2 and no outbound HTTPS, so annotation happens once on a workstation and the columns EM
    # needs are staged as compact parquet (~20 MB for all 7 loci vs 1.6 GB of raw arda TSV).
    data_dir = os.environ.get("EM_DATA_DIR")
    if data_dir:
        uniq = data.unique_clonotypes(pl.read_parquet(Path(data_dir) / f"human_{locus}_nonfunctional.parquet"))
    else:
        uniq = data.prepare("human", locus, "nonfunctional", out_dir=str(WORK))
    # Filter on GENE, not allele. arda and OLGA resolve alleles differently -- arda calls
    # TRBV20-1*07, which OLGA's 89-allele index does not contain -- and an allele-level
    # `is_in(vset)` therefore deletes the WHOLE GENE before gene_masks (below) ever sees it.
    # Measured on human TRB: TRBV20-1, the most-used human TRBV, went to 0 training clonotypes
    # and hence P(V)=0 in the shipped model. gene_masks already maps a call to all model alleles
    # of its gene, so the read is perfectly usable -- the pre-filter was throwing it away.
    # Gene-level keeps 32,562 vs 24,980 clonotypes and 54 vs 51 V genes.
    vg = pl.col("v_call").str.split("*").list.first()
    jg = pl.col("j_call").str.split("*").list.first()
    # ALL non-functional reads -- out-of-frame AND stop-codon alike. Both escaped selection, which
    # is the only property the generative model needs. Restricting to out-of-frame is itself a
    # selection: it keeps only junctions whose length happens not to be a multiple of 3, and
    # conditioning the training set on a junction-length property is exactly the kind of bias the
    # insertion-length model would then learn. Using the whole bucket is both unbiased and ~5x the
    # data (TRB: 32.6k out-of-frame vs 139.6k total unique clonotypes).
    uniq = uniq.filter(
        vg.is_in(list(vgenes)) & jg.is_in(list(jgenes))
        & pl.col("junction").str.to_uppercase().str.contains(r"^[ACGT]+$"))
    n_all = uniq.height
    # No cap and no sampling: EM runs on every clonotype that survives the germline filter. A
    # subsample is a silent statement that the tail does not matter, and the tail is where the
    # rare V genes live -- the ones that collapse to P(V)=0 in the first place.
    seqs = [s.upper() for s in uniq["junction"].to_list()]
    # SOFT REALIGNMENT (Murugan/IGoR): consider ALL V/J candidates per read and let EM's soft
    # counts determine the assignment, rather than hard-masking to arda's single best call. arda
    # commits V-choice to one best-match allele, so a functional gene it systematically under-calls
    # (TRBV12-4, TRBV11-3 -- observed solo hundreds of times, junctions scoring Pgen>0) never gets
    # soft count and is driven to P(V)=0. An empty per-read V/J mask == all candidates (verified
    # bit-identical to masks=None), so the E-step marginalises over every gene and a
    # read incompatible with a gene contributes ~0 to it via Pgen -- no gene is ever hard-zeroed.
    # ~15x slower than the arda mask (all V vs one), hence Aldan-3. D stays arda-masked: D is not
    # the zeroing problem and its short germline makes the mask a real, safe speedup.
    soft = os.environ.get("SOFT_REALIGN") == "1"
    masks = ([([], [], None) for _ in seqs] if soft          # empty V/J list == all candidates
             else gene_masks(base, uniq["v_call"].to_list(), uniq["j_call"].to_list()))
    if base.chain_type == "VDJ":
        # arda's D call, as a SET -- AIRR writes an aligner tie comma-separated
        # ("IGHD2-2*01,IGHD2-2*02,IGHD2-2*03"), so the old `r["d_call"] in dset` tested the whole
        # joined string against a set of single alleles, never matched, and fell through to an
        # empty mask = all D enumerated. Measured on IGH: 19,451/160,324 (12.1%) of D calls are
        # ambiguous, and 64.1% of reads ended up with an empty D mask. An empty mask is not just
        # slow (35 D genes for IGH vs 3 for TRB), it is a WORSE model: the D marginal stops being
        # anchored to what arda actually saw. A genuinely absent D call still yields [] -- correct,
        # since we know nothing -- and that is 50% of IGH on its own (short D + SHM).
        def d_mask(call: str | None) -> list[str]:
            if not call:
                return []
            return [a for a in (p.strip() for p in call.split(",")) if a in dset]
        masks = [(mk[0], mk[1], d_mask(r.get("d_call")))
                 for mk, r in zip(masks, uniq.iter_rows(named=True))]

    # Learn tandem-D anchored to arda: on the D-bearing loci a read may be n_D=2 only where arda
    # called a second D (d2_call). This counters the tandem-vs-long-insertion identifiability that
    # inflates unregularized D-D EM (TRB drifts to P(n_D=2)~0.28). EM_SINGLE_D=1 forces strict
    # single-D; ND_PRIOR adds a Dirichlet single-D pseudocount on top.
    single_d = os.environ.get("EM_SINGLE_D") == "1"
    nd_prior = float(os.environ.get("ND_PRIOR", "0"))
    # Dirichlet pseudocount over the germline's functional V/J alleles. P(V)=0 is an EM
    # absorbing state, so without this one unlucky iteration permanently deletes a real gene:
    # unregularized, human TRB kept 30 of OLGA's 57 V genes having SEEN 54 in the data.
    gene_prior = float(os.environ.get("GENE_PRIOR", "1.0"))
    dd_allowed = None
    if base.chain_type == "VDJ" and not single_d:
        dd_allowed = [r.get("d2_call") is not None for r in uniq.iter_rows(named=True)]
    t = time.perf_counter()
    model, rep = infer_native(base, seqs, masks=masks, max_iter=ITERS, tol=0.0,
                              single_d=single_d, dd_allowed=dd_allowed, nd_prior=nd_prior,
                              gene_prior=gene_prior)
    dt = time.perf_counter() - t
    out = DEST / locus
    model.save(out)
    nd = dict(zip(model.tables["n_d"]["n_d"].to_list(), model.tables["n_d"]["p"].to_list())) if "n_d" in model.tables else {}
    return {"locus": locus, "chain": base.chain_type, "n_clono": n_all, "n_used": len(seqs),
            "p_nd2": nd.get(2), "ll0": rep.loglik[0], "ll1": rep.loglik[-1], "iters": rep.n_iter, "sec": dt}


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    want = os.environ.get("LOCI")
    loci = {k: LOCI[k] for k in want.split(",")} if want else LOCI
    rows = []
    for locus, name in loci.items():
        print(f"[{locus}] building learned model ...", flush=True)
        try:
            r = build(locus, name)
            rows.append(r)
            print(f"  [{locus}] {r['n_used']}/{r['n_clono']} clonotypes, {r['iters']} iters, {r['sec']:.0f}s, "
                  f"P(n_D=2)={r['p_nd2']}, TRAINING LL {r['ll0']:.0f}->{r['ll1']:.0f}", flush=True)
        except Exception as e:  # noqa: BLE001 — report and continue
            print(f"  [{locus}] FAILED {type(e).__name__}: {e}", flush=True)
    print("\n| locus | chain | clonotypes | used | P(n_D=2) | TRAINING LL (0->N) | iters | sec |")
    print("|-------|-------|-----------|------|----------|--------------------|-------|-----|")
    for r in rows:
        p = f"{r['p_nd2']:.4f}" if r["p_nd2"] is not None else "—"
        print(f"| {r['locus']} | {r['chain']} | {r['n_clono']} | {r['n_used']} | {p} | "
              f"{r['ll0']:.0f} -> {r['ll1']:.0f} | {r['iters']} | {r['sec']:.0f} |")


if __name__ == "__main__":
    main()
