"""Benchmark amino-acid Pgen throughput for the vdjtools appendix (murugan_model.tex).

Compares the native ``_core`` transfer matrix (vdjtools.model.native.pgen_aa) against OLGA on the
human TRB (VDJ) and TRA (VJ) OLGA default models, single-threaded, for two query types:
  * exact  -- Pgen of one amino-acid CDR3 (native pgen_aa            vs OLGA compute_aa_CDR3_pgen),
  * 1-mm   -- Pgen of the Hamming-1 ball (native pgen_aa mismatches=1 vs OLGA compute_hamming_dist_1_pgen).
All queries marginalize over every functional V/J/D. Prints the per-query throughput table and the
LaTeX macro values the appendix substitutes into its tables. Set BENCH_ENUM=1 to also time the pure-
Python scenario enumeration on the single shortest CDR3 (it is the correctness oracle, not a runtime
path, and is O(seconds/query) unrestricted).

Run:  VDJTOOLS_OLGA_MODELS=/path/to/olga/default_models python appendix/bench_pgen.py
"""
from __future__ import annotations

import os
import time
from pathlib import Path

from vdjtools.model import from_olga, native
from vdjtools.model.pgen import pgen_aa as enum_pgen_aa
from vdjtools.model.pgen import prepare
from vdjtools.model.generate import generate

MODELS = Path(
    os.environ.get(
        "VDJTOOLS_OLGA_MODELS",
        "/Users/mikesh/vcs/code/mirpy/mir/resources/olga/default_models",
    )
)


def _olga_model(sub, vdj):
    from olga import load_model as ol
    import olga.generation_probability as gp

    d = MODELS / sub
    gen = (ol.GenerativeModelVDJ if vdj else ol.GenerativeModelVJ)()
    gen.load_and_process_igor_model(str(d / "model_marginals.txt"))
    gd = (ol.GenomicDataVDJ if vdj else ol.GenomicDataVJ)()
    gd.load_igor_genomic_data(
        str(d / "model_params.txt"),
        str(d / "V_gene_CDR3_anchors.csv"),
        str(d / "J_gene_CDR3_anchors.csv"),
    )
    return (gp.GenerationProbabilityVDJ if vdj else gp.GenerationProbabilityVJ)(gen, gd)


def _time(fn, seqs, repeats):
    t0 = time.perf_counter()
    out = []
    for _ in range(repeats):
        out = [fn(s) for s in seqs]
    return (time.perf_counter() - t0) / repeats / len(seqs) * 1e3, out  # ms / sequence


def bench(sub, locus, vdj, macros):
    m = from_olga(MODELS / sub, locus=locus)
    olga = _olga_model(sub, vdj)
    df = sorted(generate(m, 60, seed=7, productive_only=True).to_dicts(), key=lambda r: len(r["cdr3_aa"]))
    seqs = [r["cdr3_aa"] for r in df if len(r["cdr3_aa"]) <= 15][:20]
    print(f"\n=== {locus} ({'VDJ' if vdj else 'VJ'})  n={len(seqs)}  lengths={sorted({len(s) for s in seqs})} ===")

    native.pgen_aa(m, seqs[0]); native.pgen_aa(m, seqs[0], mismatches=1)  # warm pack()
    e_tm, e_out = _time(lambda s: native.pgen_aa(m, s), seqs, 5)
    e_ol, e_ref = _time(lambda s: olga.compute_aa_CDR3_pgen(s), seqs, 2)
    h_tm, h_out = _time(lambda s: native.pgen_aa(m, s, mismatches=1), seqs, 3)
    h_ol, h_ref = _time(lambda s: olga.compute_hamming_dist_1_pgen(s), seqs, 1)
    e_rel = max(abs(a - b) / b for a, b in zip(e_out, e_ref) if b > 0)
    h_rel = max(abs(a - b) / b for a, b in zip(h_out, h_ref) if b > 0)

    print(f"{'query':<7}{'OLGA ms':>10}{'native ms':>12}{'speedup':>10}{'max rel':>12}")
    print(f"{'exact':<7}{e_ol:>10.3f}{e_tm:>12.3f}{e_ol / e_tm:>9.1f}x{e_rel:>12.1e}")
    print(f"{'1-mm':<7}{h_ol:>10.3f}{h_tm:>12.3f}{h_ol / h_tm:>9.1f}x{h_rel:>12.1e}")
    macros[locus] = dict(e_ol=e_ol, e_tm=e_tm, e_x=e_ol / e_tm, h_ol=h_ol, h_tm=h_tm, h_x=h_ol / h_tm)

    if os.getenv("BENCH_ENUM") and vdj:
        prep = prepare(m)
        enum_ms, _ = _time(lambda s: enum_pgen_aa(prep, s), [seqs[0]], 1)
        print(f"{'enum':<7}{'':>10}{enum_ms:>12.0f}  (single L={len(seqs[0])} exact query — oracle only)")
        macros["enum"] = enum_ms


def main():
    macros = {}
    bench("human_T_beta", "TRB", True, macros)
    bench("human_T_alpha", "TRA", False, macros)
    b = macros["TRB"]
    print("\n--- LaTeX macro values ---")
    print(f"TRB exact: OLGA {b['e_ol']:.2f}  native {b['e_tm']:.2f}  {b['e_x']:.1f}x")
    print(f"TRB 1-mm : OLGA {b['h_ol']:.1f}  native {b['h_tm']:.1f}  {b['h_x']:.1f}x")
    a = macros["TRA"]
    print(f"TRA 1-mm : OLGA {a['h_ol']:.1f}  native {a['h_tm']:.1f}  {a['h_x']:.1f}x")


if __name__ == "__main__":
    main()
