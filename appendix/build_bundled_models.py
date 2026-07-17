"""Build the EM-learned bundled models from real HuggingFace out-of-frame reads, all 7 loci.  2026-07-12

For each locus: fetch nonfunctional reads (HF), arda-map to unique clonotypes, and run native EM
(D-D by default on the D-bearing loci IGH/TRD/TRB) seeded from the OLGA bootstrap model. Saves to
python/vdjtools/model/_bundled/learned/<locus>/ (parquet marginals + manifest.json).

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
OLGA = Path(os.environ.get("VDJTOOLS_OLGA_MODELS", str(Path(_olga.__file__).parent / "default_models")))
DEST = Path("python/vdjtools/model/_bundled/learned")
# NB the previous default pointed into ANOTHER SESSION's scratchpad and the build silently
# reused a truncated arda cache from it -- that is how the shipped TRB model came to be
# trained on 870 clonotypes when 32,562 were available. arda.annotate_reads does its own
# caching under out_dir; there is no second cache layer here any more.
WORK = Path(os.environ.get("EM_WORK", "/tmp/em_work"))
LOCI = {"TRA": "human_T_alpha", "TRB": "human_T_beta", "TRG": "human_T_gamma", "TRD": "human_T_delta",
        "IGH": "human_B_heavy", "IGK": "human_B_kappa", "IGL": "human_B_lambda"}
CAP = int(os.environ.get("EM_CAP", "20000"))
ITERS = int(os.environ.get("EM_ITERS", "12"))


def build(locus: str, name: str) -> dict:
    base = from_olga(OLGA / name, locus=locus)
    dset = set(base.genomic["genes_d"]["d_allele"].to_list()) if base.chain_type == "VDJ" else set()
    vgenes = set(_gene_to_alleles(base, "v"))
    jgenes = set(_gene_to_alleles(base, "j"))

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
    uniq = uniq.filter(
        vg.is_in(list(vgenes)) & jg.is_in(list(jgenes))
        & pl.col("junction").str.to_uppercase().str.contains(r"^[ACGT]+$")
        # OUT-OF-FRAME ONLY, by junction length -- the Murugan/OLGA training convention. The HF
        # `nonfunctional` bucket is out-of-frame + stop-codon, and the stop half is 68% of it and
        # NOT equivalent: 56.5% of those reads carry the stop INSIDE the junction, so training on
        # them conditions the junction's nucleotide composition on containing a stop -- exactly
        # what the insertion/dinucleotide model learns. Length mod 3 != 0 is used rather than
        # arda's vj_in_frame flag because the two disagree on ~16k TRB reads, and "never
        # translated" is a property of the length, not of a flag.
        & ((pl.col("junction").str.len_bytes() % 3) != 0))
    n_all = uniq.height
    if n_all > CAP:
        uniq = uniq.sample(CAP, seed=1)
    seqs = [s.upper() for s in uniq["junction"].to_list()]
    masks = gene_masks(base, uniq["v_call"].to_list(), uniq["j_call"].to_list())
    if base.chain_type == "VDJ":  # add the arda-aligned D to each read's mask
        masks = [(mk[0], mk[1], [r["d_call"]] if r.get("d_call") in dset else [])
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
