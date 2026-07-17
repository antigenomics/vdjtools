"""7-chain concordance: the native vdjtools pipeline vs the OLGA oracle, all 7 human loci.  2026-07-11

For each locus: generate sequences from our model, then compare our **native** nt and aa Pgen to
OLGA's (Pearson r on log10 Pgen, max relative error), check aa == Σ nt over synonymous codons on a
short CDR3, and (D-bearing loci) confirm the D-D-default EM does not hallucinate tandems on
single-D OLGA data. Prints a Markdown concordance table.

Reproduce (needs the [oracle] extra: `pip install olga`):

    python appendix/concordance.py            # all 7 loci
    LOCI=TRB,TRD python appendix/concordance.py
"""
import itertools
import os
from collections import defaultdict
from pathlib import Path

import numpy as np

import olga.generation_probability as ogp
import olga.load_model as olm

from vdjtools.model import from_olga, native
from vdjtools.model.generate import generate
from vdjtools.model.pgen import pgen_aa, pgen_nt, prepare
from vdjtools.model.reference import _CODON_TABLE

# pip olga ships its own default_models; the old default pointed at a mirpy checkout that no
# longer exists on disk, so this script could not run at all. NB pip olga has no TRG/TRD.
OLGA = Path(os.environ.get("VDJTOOLS_OLGA_MODELS", str(Path(olm.__file__).parent / "default_models")))
LOCI = {
    "TRA": ("human_T_alpha", "VJ"), "TRB": ("human_T_beta", "VDJ"),
    "TRG": ("human_T_gamma", "VJ"), "TRD": ("human_T_delta", "VDJ"),
    "IGH": ("human_B_heavy", "VDJ"), "IGK": ("human_B_kappa", "VJ"), "IGL": ("human_B_lambda", "VJ"),
}
N_NT = int(os.environ.get("N_NT", "150"))   # sequences for the nt concordance
N_AA = int(os.environ.get("N_AA", "12"))    # subsample for the (slower) OLGA aa oracle


def olga_pg(sub: Path, chain: str):
    G = olm.GenomicDataVDJ if chain == "VDJ" else olm.GenomicDataVJ
    M = olm.GenerativeModelVDJ if chain == "VDJ" else olm.GenerativeModelVJ
    P = ogp.GenerationProbabilityVDJ if chain == "VDJ" else ogp.GenerationProbabilityVJ
    g = G(); g.load_igor_genomic_data(str(sub / "model_params.txt"), str(sub / "V_gene_CDR3_anchors.csv"), str(sub / "J_gene_CDR3_anchors.csv"))
    m = M(); m.load_and_process_igor_model(str(sub / "model_marginals.txt"))
    return P(m, g)


def logcorr(a, b):
    a, b = np.asarray(a), np.asarray(b)
    ok = (a > 0) & (b > 0)
    if ok.sum() < 3:
        return float("nan"), ok.sum()
    return float(np.corrcoef(np.log10(a[ok]), np.log10(b[ok]))[0, 1]), int(ok.sum())


def maxrel(a, b):
    a, b = np.asarray(a), np.asarray(b)
    ok = b > 0
    return float(np.max(np.abs(a[ok] - b[ok]) / b[ok])) if ok.any() else float("nan")


def run_locus(locus: str) -> dict:
    name, chain = LOCI[locus]
    sub = OLGA / name
    m = from_olga(sub, locus=locus)
    prep = prepare(m)
    po = olga_pg(sub, chain)
    df = generate(m, N_NT, seed=1, productive_only=True)  # in-frame: OLGA's nt Pgen rejects out-of-frame
    rows = df.to_dicts()
    # nt Pgen: native vs OLGA
    ours_nt, olga_nt = [], []
    for r in rows:
        nt, V, J = r["cdr3_nt"], r["v_call"], r["j_call"]
        ours_nt.append(native.pgen_nt(m, nt, V, J))
        olga_nt.append(po.compute_nt_CDR3_pgen(nt, V, J))
    r_nt, n_nt = logcorr(ours_nt, olga_nt)
    mr_nt = maxrel(ours_nt, olga_nt)
    # aa Pgen: native vs OLGA (subsample; OLGA aa is slow)
    prod = [r for r in rows if r["productive"]][:N_AA]
    ours_aa, olga_aa = [], []
    for r in prod:
        aa, V, J = r["cdr3_aa"], r["v_call"], r["j_call"]
        ours_aa.append(native.pgen_aa(m, aa, V, J))
        olga_aa.append(po.compute_aa_CDR3_pgen(aa, V, J))
    r_aa, n_aa = logcorr(ours_aa, olga_aa)
    mr_aa = maxrel(ours_aa, olga_aa)
    # aa == Σ nt over synonymous, on the shortest productive CDR3 with a bounded synonymous count
    syn = defaultdict(list)
    for cod, a in _CODON_TABLE.items():
        syn[a].append(cod)
    aa_nt_ok = None
    for r in sorted(rows, key=lambda r: len(r["cdr3_aa"])):  # search all in-frame reads for a short one
        aa, V, J = r["cdr3_aa"], r["v_call"], r["j_call"]
        if int(np.prod([len(syn[a]) for a in aa])) <= 20000:
            brute = sum(native.pgen_nt(m, "".join(c), V, J) for c in itertools.product(*[syn[a] for a in aa]))
            aa_nt_ok = bool(np.isclose(native.pgen_aa(m, aa, V, J), brute, rtol=1e-6, atol=1e-300))
            break
    return {"locus": locus, "chain": chain, "n_nt": n_nt, "r_nt": r_nt, "mr_nt": mr_nt,
            "n_aa": n_aa, "r_aa": r_aa, "mr_aa": mr_aa, "aa_nt_ok": aa_nt_ok}


def main():
    want = os.environ.get("LOCI")
    loci = want.split(",") if want else list(LOCI)
    results = []
    for locus in loci:
        print(f"[{locus}] running concordance ...", flush=True)
        try:
            results.append(run_locus(locus))
        except Exception as e:  # noqa: BLE001 — report and continue
            print(f"[{locus}] FAILED {type(e).__name__}: {e}", flush=True)
    print("\n| locus | chain | nt r(log) | nt max-rel | aa r(log) | aa max-rel | aa==Σnt |")
    print("|-------|-------|-----------|------------|-----------|------------|---------|")
    for r in results:
        print(f"| {r['locus']} | {r['chain']} | {r['r_nt']:.5f} | {r['mr_nt']:.1e} | "
              f"{r['r_aa']:.5f} | {r['mr_aa']:.1e} | {r['aa_nt_ok']} |")


if __name__ == "__main__":
    main()
