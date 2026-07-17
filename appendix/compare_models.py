"""Three-way model comparison: OLGA (pip) vs OLGA-in-vdjtools vs the EM-learned model.  2026-07-17

Three questions, deliberately kept apart -- conflating them is how a broken model ships next to a
correct one and looks like its peer.

1. CORRECTNESS (must be exact). `vdjtools`'s native Pgen on the bundled `olga` model must equal
   the `olga` pip package's own Pgen to machine precision, on the same model files. This is an
   invariant, not a metric: any deviation is a bug in our DP. OLGA is the oracle here.

2. THE JUNCTION MODEL (the real comparison). Trims, insertion lengths, insertion dinucleotides,
   D usage, n_D -- the events that describe how a junction is BUILT. These are what the HF
   out-of-frame reads are good for, and where `learned` and `olga` can be meaningfully compared:
   both are estimates of the same physical recombination machinery, so large divergence is a
   finding either way.

3. V/J USAGE (NOT a defect when it differs). P(V) and P(J|V) are PROTOCOL-DEPENDENT: the HF
   reads are 5'RACE, OLGA's model was fit to DNA multiplex (Adaptive). The two amplify different
   V's at different rates, so the usage marginals SHOULD disagree and neither is "wrong" -- e.g.
   TRBV19 is 37% of these functional reads vs OLGA's 3.1%. The usable content here is that every
   V/J combination is represented, so J|V and the junction model are learnable; the usage itself
   is meant to be rescaled to whatever protocol the user's own sample came from
   (`vdjtools.model.rescale_usage`). What IS a defect is a real gene at P(V)=0 -- that makes Pgen
   exactly 0 for every clonotype using it, silently, and no correlation summary would show it.
   Reported per GENE, not per allele: the arda allele split (TRBV19*03 vs *01) is mismapping-level
   noise on 151bp reads, so allele-resolution usage is not a meaningful quantity from this data.

Reproduce (needs the [oracle] extra + HF access to isalgo/airr_model_read):

    python appendix/compare_models.py                 # TRB
    LOCI=TRB,TRA python appendix/compare_models.py
"""
import os
from pathlib import Path

import numpy as np
import olga as _olga
import olga.generation_probability as ogp
import olga.load_model as olm
import polars as pl

from vdjtools.model import Model, data, native
from vdjtools.model.schema import normalization_keys

# Load the models from THIS REPO's tree, not via bundled.load_bundled. load_bundled resolves
# through the installed package, which under an editable install can point at a different
# checkout entirely -- so a freshly rebuilt model would be written here and silently NOT read
# back, and the comparison would report the old model's numbers as the new one's. The builder
# writes to this same relative path; keep the two in agreement.
BUNDLED = Path("python/vdjtools/model/_bundled")
# The OLGA models shipped in THIS repo, not pip olga's: pip ships only 5 human loci (no
# TRG/TRD) plus mouse, while tests/python/fixtures/olga/default_models carries all 7 human loci.
# The TRG/TRD marginals originate from mirpy's legacy-v2 branch (commit aeccd75) and are verified
# byte-identical to what the bundled parquet were built from; olga-pip scores with them fine, so
# they are a real oracle for those two loci, which pip alone cannot be.
_REPO_OLGA = Path(__file__).resolve().parent.parent / "tests" / "python" / "fixtures" / "olga" / "default_models"
OLGA = Path(os.environ.get("VDJTOOLS_OLGA_MODELS", str(_REPO_OLGA)))
WORK = Path(os.environ.get("EM_WORK", "/tmp/em_work"))
LOCI = {"TRA": ("human_T_alpha", "VJ"), "TRB": ("human_T_beta", "VDJ"),
        "IGH": ("human_B_heavy", "VDJ"), "IGK": ("human_B_kappa", "VJ"),
        "IGL": ("human_B_lambda", "VJ")}          # pip olga ships no TRG/TRD
N = int(os.environ.get("N_CMP", "300"))


def olga_pgen_model(sub: Path, chain: str):
    G = olm.GenomicDataVDJ if chain == "VDJ" else olm.GenomicDataVJ
    M = olm.GenerativeModelVDJ if chain == "VDJ" else olm.GenerativeModelVJ
    P = ogp.GenerationProbabilityVDJ if chain == "VDJ" else ogp.GenerationProbabilityVJ
    g = G()
    g.load_igor_genomic_data(str(sub / "model_params.txt"), str(sub / "V_gene_CDR3_anchors.csv"),
                             str(sub / "J_gene_CDR3_anchors.csv"))
    m = M()
    m.load_and_process_igor_model(str(sub / "model_marginals.txt"))
    return P(m, g)


def _per_parent_tv(a: pl.DataFrame, b: pl.DataFrame, parents: list[str]) -> list[float]:
    """TV distance between two marginal tables, one value per conditioning-value group.

    Each group is a proper probability distribution summing to 1 (or to 0 where the model has no
    mass for that parent), so its TV is in [0,1] and is interpretable on its own. Groups where
    EITHER model has zero mass are skipped rather than scored 1: a V that OLGA never uses says
    nothing about how ITS junctions are trimmed, it only says the two protocols use different V's.
    """
    value_cols = [c for c in a.columns if c != "p" and c not in parents]
    if not parents:                         # unconditional table (vd_ins, j_choice…) -> one group
        da = {tuple(r[:-1]): r[-1] for r in a.select([*value_cols, "p"]).iter_rows()}
        db = {tuple(r[:-1]): r[-1] for r in b.select([*value_cols, "p"]).iter_rows()}
        return [0.5 * sum(abs(da.get(k, 0.0) - db.get(k, 0.0)) for k in set(da) | set(db))]
    out = []
    ga = {k: g for k, g in a.group_by(parents)}
    gb = {k: g for k, g in b.group_by(parents)}
    for key in set(ga) & set(gb):
        fa, fb = ga[key], gb[key]
        if fa["p"].sum() <= 0 or fb["p"].sum() <= 0:
            continue                        # a parent one model never uses — see docstring
        da = {tuple(r[:-1]): r[-1] for r in fa.select([*value_cols, "p"]).iter_rows()}
        db = {tuple(r[:-1]): r[-1] for r in fb.select([*value_cols, "p"]).iter_rows()}
        out.append(0.5 * sum(abs(da.get(k, 0.0) - db.get(k, 0.0)) for k in set(da) | set(db)))
    return out


def _entropy(p: np.ndarray) -> float:
    """Shannon entropy in bits of a (possibly unnormalized/multi-group) probability column."""
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    p = p / p.sum()
    return float(-np.sum(p * np.log2(p)))


def v_support(model) -> dict:
    """Which V GENES carry non-zero generative mass, and how concentrated is it?

    Per gene, not per allele: on 151bp 5'RACE reads arda's allele split (TRBV19*03 vs *01) is
    mismapping-level noise, so allele-resolution usage is not a meaningful quantity here.
    """
    t = model.tables["v_choice"]
    g = (t.with_columns(pl.col("v_allele").str.split("*").list.first().alias("g"))
          .group_by("g").agg(pl.col("p").sum().alias("p")))
    p = g["p"].to_numpy()
    top = g.sort("p", descending=True).head(1)
    return {"n_g": g.height, "n_g_nonzero": int((p > 0).sum()),
            "top_v": top["g"][0], "top_p": float(top["p"][0])}


def held_out_loglik(model, seqs: list[str], v_calls: list[str], j_calls: list[str]) -> tuple:
    """Mean log10 Pgen over held-out out-of-frame junctions, and how many score exactly 0.

    A Pgen of 0 is not a small number -- it is the model asserting the sequence is impossible.
    Those are excluded from the mean (log10(0) = -inf) and counted separately, because that count
    IS the finding when a model has collapsed V support.
    """
    from vdjtools.model.infer import _gene_to_alleles
    va, ja = _gene_to_alleles(model, "v"), _gene_to_alleles(model, "j")
    out = []
    for s, v, j in zip(seqs, v_calls, j_calls):
        vs, js = va.get(v.split("*")[0], []), ja.get(j.split("*")[0], [])
        if not vs or not js:
            out.append(0.0)
            continue
        tot = 0.0
        for vv in vs:
            for jj in js:
                tot += native.pgen_nt(model, s.upper(), vv, jj)
        out.append(tot)
    a = np.array(out)
    nz = a > 0
    return float(np.mean(np.log10(a[nz]))) if nz.any() else float("nan"), int((~nz).sum()), len(a)


def run(locus: str, name: str, chain: str) -> None:
    print(f"\n{'='*78}\n{locus}\n{'='*78}", flush=True)
    m_olga = Model.load(BUNDLED / "olga" / locus)
    m_learn = Model.load(BUNDLED / "learned" / locus)
    print(f"  olga    <- {BUNDLED / 'olga' / locus}")
    print(f"  learned <- {BUNDLED / 'learned' / locus}")
    opg = olga_pgen_model(OLGA / name, chain)

    # ---- test sequences: real held-out FUNCTIONAL junctions (what a user actually scores) ----
    fq_reads = WORK / f"{locus}_functional_reads.parquet"
    if fq_reads.exists():
        reads = pl.read_parquet(fq_reads)
    else:
        reads = data.annotate_reads(data.fetch_fastq("human", locus, "functional"), out_dir=str(WORK),
                                    prefix=f"human_{locus}_functional", organism="human", cap=20000)
        reads.write_parquet(fq_reads)
    uniq = data.unique_clonotypes(reads).filter(
        pl.col("junction").str.to_uppercase().str.contains(r"^[ACGT]+$")
        & ((pl.col("junction").str.len_bytes() % 3) == 0))
    test = uniq.head(N)
    seqs = [s.upper() for s in test["junction"].to_list()]

    # ---- 1. CORRECTNESS: vdjtools-on-olga-model MUST equal olga-pip ----
    ours = np.array([native.pgen_nt(m_olga, s, None, None) for s in seqs])
    theirs = np.array([opg.compute_nt_CDR3_pgen(s) for s in seqs])
    both = (ours > 0) & (theirs > 0)
    rel = np.abs(ours[both] - theirs[both]) / theirs[both]
    r = np.corrcoef(np.log10(ours[both]), np.log10(theirs[both]))[0, 1]
    print(f"[1] CORRECTNESS  vdjtools(olga model) vs olga-pip, {both.sum()}/{len(seqs)} both>0")
    print(f"    r(log10 Pgen) = {r:.10f}    max rel err = {rel.max():.3e}"
          f"    -> {'EXACT' if rel.max() < 1e-9 else 'MISMATCH — BUG'}")

    # ---- 2. THE JUNCTION MODEL: where learned and olga are genuinely comparable ----
    print(f"\n[2] JUNCTION MODEL — mean per-parent total-variation distance, learned vs olga")
    print(f"    (these describe how a junction is BUILT; both models estimate the same machinery)")
    print(f"    TV is in [0,1]: 0 = identical, 1 = disjoint. Conditional events are averaged over")
    print(f"    their conditioning value, and ONLY over parents where BOTH models have mass --")
    print(f"    a parent OLGA never uses is a usage difference (section 3), not a junction one.")
    print(f"    {'event':14s} {'mean TV':>8s} {'n parents':>10s} {'H(olga)':>9s} {'H(learned)':>11s}")
    for ev in ("v_3_del", "j_5_del", "d_del", "vd_ins", "dj_ins", "vd_dinucl", "dj_dinucl",
               "d_gene", "j_choice", "vj_ins", "vj_dinucl"):
        if ev not in m_olga.tables or ev not in m_learn.tables:
            continue
        a, b = m_olga.tables[ev], m_learn.tables[ev]
        event = m_olga.manifest.events[ev]
        parents = list(normalization_keys(event))
        tvs = _per_parent_tv(a, b, parents)
        if not tvs:
            continue
        ha, hb = _entropy(a["p"].to_numpy()), _entropy(b["p"].to_numpy())
        print(f"    {ev:14s} {np.mean(tvs):>8.4f} {len(tvs):>10d} {ha:>9.3f} {hb:>11.3f}")

    # ---- 3. V/J USAGE: protocol-dependent, reported per GENE ----
    print(f"\n[3] V-GENE USAGE — protocol-dependent (5'RACE here vs DNA-multiplex for OLGA).")
    print(f"    Divergence is EXPECTED and is rescaled per-sample; P(V)=0 on a real gene is NOT.")
    print(f"    {'model':10s} {'V genes':>9s} {'nonzero':>8s} {'top gene':>12s} {'top mass':>9s}")
    for lbl, mm in (("olga", m_olga), ("learned", m_learn)):
        s = v_support(mm)
        print(f"    {lbl:10s} {s['n_g']:>9d} {s['n_g_nonzero']:>8d} {s['top_v']:>12s} {s['top_p']:>9.4f}")

    oof_p = WORK / f"{locus}_oof_heldout.parquet"
    if oof_p.exists():
        oof = pl.read_parquet(oof_p)
        print(f"\n[4] HELD-OUT LOG-LIKELIHOOD on {oof.height:,} out-of-frame junctions the EM never saw")
        print(f"    NB this is confounded by (3): the learned model has this protocol's V usage and")
        print(f"    the held-out set is from the same protocol, so it is favoured on usage alone.")
        print(f"    {'model':10s} {'mean log10 Pgen':>16s} {'Pgen==0':>9s} {'n':>7s}")
        for lbl, mm in (("olga", m_olga), ("learned", m_learn)):
            ll, nz, n = held_out_loglik(mm, oof["junction"].to_list(), oof["v_call"].to_list(),
                                        oof["j_call"].to_list())
            print(f"    {lbl:10s} {ll:>16.4f} {nz:>9d} {n:>7d}")


def main() -> None:
    want = os.environ.get("LOCI", "TRB")
    for locus in want.split(","):
        name, chain = LOCI[locus]
        run(locus, name, chain)


if __name__ == "__main__":
    main()
