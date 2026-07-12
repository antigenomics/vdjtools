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
from vdjtools.model.infer import gene_masks, infer_native

OLGA = Path(os.environ.get("VDJTOOLS_OLGA_MODELS", "/Users/mikesh/vcs/code/mirpy/mir/resources/olga/default_models"))
DEST = Path("python/vdjtools/model/_bundled/learned")
WORK = Path(os.environ.get("EM_WORK", "/private/tmp/claude-501/-Users-mikesh-vcs-code-vdjtools/84fb4c31-220b-41cf-a8f1-e4fce9fefd56/scratchpad/_bundled_reads"))
LOCI = {"TRA": "human_T_alpha", "TRB": "human_T_beta", "TRG": "human_T_gamma", "TRD": "human_T_delta",
        "IGH": "human_B_heavy", "IGK": "human_B_kappa", "IGL": "human_B_lambda"}
CAP = int(os.environ.get("EM_CAP", "2000"))
ITERS = int(os.environ.get("EM_ITERS", "12"))


def build(locus: str, name: str) -> dict:
    base = from_olga(OLGA / name, locus=locus)
    vset = set(base.genomic["genes_v"][f"{'v'}_allele"])
    jset = set(base.genomic["genes_j"]["j_allele"])
    dset = set(base.genomic["genes_d"]["d_allele"].to_list()) if base.chain_type == "VDJ" else set()

    airr = WORK / f"human_{locus}_nonfunctional.airr.tsv"  # reuse arda's output if already mapped
    if airr.exists():
        uniq = data.unique_clonotypes(pl.read_csv(airr, separator="\t", infer_schema_length=20000))
    else:
        uniq = data.prepare("human", locus, "nonfunctional", out_dir=str(WORK))
    uniq = uniq.filter(
        pl.col("v_call").is_in(list(vset)) & pl.col("j_call").is_in(list(jset))
        & pl.col("junction").str.to_uppercase().str.contains(r"^[ACGT]+$"))
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
    dd_allowed = None
    if base.chain_type == "VDJ" and not single_d:
        dd_allowed = [r.get("d2_call") is not None for r in uniq.iter_rows(named=True)]
    t = time.perf_counter()
    model, rep = infer_native(base, seqs, masks=masks, max_iter=ITERS, tol=0.0,
                              single_d=single_d, dd_allowed=dd_allowed, nd_prior=nd_prior)
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
                  f"P(n_D=2)={r['p_nd2']}, held-out LL {r['ll0']:.0f}->{r['ll1']:.0f}", flush=True)
        except Exception as e:  # noqa: BLE001 — report and continue
            print(f"  [{locus}] FAILED {type(e).__name__}: {e}", flush=True)
    print("\n| locus | chain | clonotypes | used | P(n_D=2) | held-out LL (0->N) | iters | sec |")
    print("|-------|-------|-----------|------|----------|--------------------|-------|-----|")
    for r in rows:
        p = f"{r['p_nd2']:.4f}" if r["p_nd2"] is not None else "—"
        print(f"| {r['locus']} | {r['chain']} | {r['n_clono']} | {r['n_used']} | {p} | "
              f"{r['ll0']:.0f} -> {r['ll1']:.0f} | {r['iters']} | {r['sec']:.0f} |")


if __name__ == "__main__":
    main()
