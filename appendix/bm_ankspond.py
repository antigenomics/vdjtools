"""The AS/B27 motif: is it disease, or is it carriage? And can we find it without being told?

Komech 2018 (Rheumatology 57:1097) found a public TRBV9/TRBJ2-3 CDR3b motif in ankylosing
spondylitis. They found it by extending the Murugan 2012 generative model with MONTE-CARLO --
2e9 TRBV9-TRBJ2-3 rearrangements drawn to *estimate* P(amino acid sequence), then a normal fit
to log10 Pgen per sharing bin, threshold P<0.001. vdjtools computes that same quantity EXACTLY
by DP (appendix/murugan_model.tex), so every number here is the exact instrument, not a sampled
one -- and the 1-mismatch ball their MC could not reach falls out of the same call.

The confound that defines this dataset: B27 is 26/27 confounded with AS (only donor `Kal` is
AS/B27-). So "AS vs healthy" cannot separate disease from carriage. The only contrast that can
is B27-MATCHED: AS/B27+ vs healthy/B27+. That is the paper's actual claim (7 of 8 clonotypes
absent from B27+ *healthy* blood), and it is the primary here.

A second confound the metadata hides: `proj`. Five of the twelve HD/B27+ are foreign batches
(G-CSF-stimulated, CD45RA-depleted, pediatric) while every AS donor is proj=nan. The
batch-matched contrast is the honest headline; the pooled one is reported beside it.

Unit is the DONOR, never the clonotype and never the read -- Emerson tested abundance-weighting
head-to-head and it lost; weighting a 2x2 by reads is pseudoreplication (Hurlbert 1984). Four
donors (Dv/Kal/Mikh/Shep) have three samples each (PB_F + PB_CD8 + SFCD8); they collapse to one
row.

Arms:
  cohort  -- load + the 2x2 + depth check (is the contrast depth-confounded? no.)
  fig1    -- reproduce Komech's Fig.1 table from old/ (the 2018 cohort, via the new read_mitcr)
  confirm -- A: pre-registered 9-clonotype hypothesis. One test, no multiple-testing burden.
  pgen    -- B: Komech's probabilistic screen, exact DP instead of 2e9 MC draws.
  screen  -- C: V/J-pinned discovery -- can we find it WITHOUT being told the sequences?
  oracle  -- VDJdb cross-check (Yang 2022 Nature; NB partly circular -- same group/donors).

2026-07-17.
"""
from __future__ import annotations

import argparse
import resource
import time
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import fisher_exact, mannwhitneyu

from vdjtools import io as vio
from vdjtools.biomarker import stats as bstats

# --- the ground truth ---------------------------------------------------------------------------
# Komech 2018 prints NO regex and NO consensus string. The motif IS these enumerated clonotypes:
# eight from the Fig.1 table (p.1100) plus CASSAGLYSTDTQYF, named in prose on p.1101 (shared by
# all three CD8+ SF samples of B27+ patients). Cite the list, not a pattern.
MOTIF_8 = [  # Fig.1, with the paper's own (#AS, #healthy) donor counts
    ("CASSVGLYSTDTQYF", 12, 0), ("CASSVGLFSTDTQYF", 7, 0), ("CASSVGVYSTDTQYF", 6, 0),
    ("CASSVATYSTDTQYF", 5, 0), ("CASSLGLFSTDTQYF", 4, 0), ("CASSAGLFSTDTQYF", 4, 0),
    ("CASSPGLFSTDTQYF", 4, 0), ("CASSVGGFGDTQYF", 3, 1),
]
MOTIF_SF = "CASSAGLYSTDTQYF"          # p.1101, prose
MOTIF_VDJDB = "CASSVGTYSTDTQYF"       # a 10th member: in VDJdb (Yang 2022), NOT in Komech's 9.
MOTIF = [c for c, _, _ in MOTIF_8] + [MOTIF_SF]
# Our derived consensus -- NOT the authors'. Covers the eight 15-mers; the 14-mer CASSVGGFGDTQYF
# is length-incompatible and must stay enumerated. mirpy's CASS[A-Z]G[LV][YF]STDTQYF misses
# CASSVATYSTDTQYF (needs G at 6, has A) and the 14-mer -- 2 of 9.
MOTIF_RE = r"CASS[VLAP][GA][LVT][YF]STDTQYF"
V_GENE, J_GENE = "TRBV9", "TRBJ2-3"
V_ALLELE, J_ALLELE = "TRBV9*01", "TRBJ2-3*01"   # the model is keyed by ALLELE (see native._gene_idx)

REPO = "isalgo/airr_ankspond"
DEF_VDJDB = "/Users/mikesh/vcs/code/vdjdb-iedb-concordance/vdjdb_dump_2026/vdjdb.slim.txt"
RES = Path(__file__).parent / "results"


def _rss_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024**3)  # macOS: bytes


class Timer:
    def __init__(self, label): self.label = label
    def __enter__(self): self.t = time.perf_counter(); return self
    def __exit__(self, *a):
        print(f"  [{self.label}] {time.perf_counter()-self.t:.1f}s  peak RSS {_rss_gb():.2f} GB",
              flush=True)


def _snapshot() -> Path:
    from huggingface_hub import snapshot_download
    return Path(snapshot_download(REPO, repo_type="dataset", max_workers=8))


# --- cohort -------------------------------------------------------------------------------------

def load_new(snap: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    """`new/` (Komech 2022) -> (per-sample clonotypes tagged with donor_id, per-donor metadata).

    `nan` is the literal missing marker. One `Koc.tsv.gz` carries Float64 counts where every other
    file is Int64, so a naive read_csv+concat raises -- `vio.read` normalizes it away.
    """
    meta = pl.read_csv(snap / "new" / "metadata.tsv", separator="\t", infer_schema_length=0,
                       null_values=["nan"])
    frames = []
    for r in meta.iter_rows(named=True):
        p = snap / "new" / f"{r['sample_name']}.tsv.gz"
        frames.append(vio.read(str(p)).with_columns(
            pl.lit(r["donor_id"]).alias("donor"),
            pl.lit(r["sample_name"]).alias("sample"),
            pl.lit(r["fraction"]).alias("fraction"),
            pl.lit(r["sample_type"]).alias("sample_type"),
        ))
    clones = pl.concat(frames, how="vertical_relaxed")
    # donor-level design: one row per donor. `proj` is the batch column.
    donors = (meta.select("donor_id", "disease_status", "b27", "proj").unique("donor_id")
              .rename({"donor_id": "donor"}).sort("donor"))
    return clones, donors


def carriage(clones: pl.DataFrame, seqs: list[str], *, pin_v: bool = True) -> set[str]:
    """Donors carrying >=1 of `seqs`. V-pinning is load-bearing: unpinned, wrong-V convergents
    (e.g. CASSLGLFSTDTQYF on TRBV5-1/5-4/5-5/7-7/7-8) leak in and inflate the healthy arm."""
    f = pl.col("junction_aa").is_in(seqs)
    if pin_v:
        f = f & (pl.col("v_call") == V_GENE) & (pl.col("j_call") == J_GENE)
    return set(clones.filter(f)["donor"].unique().to_list())


def _fisher(carriers: set[str], pos: list[str], neg: list[str]) -> tuple[int, int, float, float]:
    a = len(carriers & set(pos)); b = len(carriers & set(neg))
    c = len(pos) - a; d = len(neg) - b
    orr, p = fisher_exact([[a, c], [b, d]], alternative="greater")
    return a, b, orr, p


def arm_cohort(clones, donors) -> None:
    print("\n=== the 2x2: B27 is 26/27 confounded with AS ===", flush=True)
    print(donors.group_by("disease_status", "b27").len().sort("disease_status", "b27"), flush=True)
    print("\n  => AS-vs-healthy cannot separate disease from carriage. B27-matched is the only "
          "contrast that can.", flush=True)
    print("\n=== `proj` batch confound among the B27+ controls ===", flush=True)
    hd = donors.filter((pl.col("disease_status") == "hd") & (pl.col("b27") == "pos"))
    print(hd.group_by("proj").len().sort("len", descending=True), flush=True)
    print(f"  AS donors with a non-null proj: "
          f"{donors.filter((pl.col('disease_status')=='as') & pl.col('proj').is_not_null()).height}"
          f"  => foreign-batch controls are not exchangeable with the cases", flush=True)

    depth = (clones.group_by("donor").agg(pl.len().alias("clones"))
             .join(donors, on="donor", how="inner"))
    print("\n=== depth: is the primary contrast depth-confounded? ===", flush=True)
    a = depth.filter((pl.col("disease_status") == "as") & (pl.col("b27") == "pos"))["clones"]
    b = depth.filter((pl.col("disease_status") == "hd") & (pl.col("b27") == "pos"))["clones"]
    u, p = mannwhitneyu(a.to_numpy(), b.to_numpy(), alternative="two-sided")
    print(f"  median clones/donor  AS/B27+ {a.median():,.0f}   HD/B27+ {b.median():,.0f}"
          f"   Mann-Whitney p={p:.2f}", flush=True)
    print("  => balanced; do NOT downsample the primary (it would only shed power).", flush=True)


def arm_fig1(snap: Path) -> None:
    """Reproduce Komech Fig.1 from `old/` -- the ACTUAL 2018 cohort, readable only via read_mitcr."""
    print("\n=== Komech 2018 Fig.1, recomputed from old/ (needs the new `mitcr` reader) ===",
          flush=True)
    meta = pl.read_csv(snap / "old" / "metadata.tsv", separator="\t", infer_schema_length=0)
    rows = []
    for r in meta.iter_rows(named=True):
        df = vio.read(str(snap / "old" / f"{r['id']}.txt.gz"))
        for m in df.filter(pl.col("junction_aa").is_in(MOTIF)).iter_rows(named=True):
            rows.append({"id": r["id"], "state": r["state"], "b27": r["b27"],
                         "cdr3": m["junction_aa"], "v": m["v_call"]})
    h = pl.DataFrame(rows)
    print(f"  {'CDR3':<17}{'#AS':>5}{'#h':>4}   {'paper':>7}", flush=True)
    exact = 0
    for c, pas, ph in MOTIF_8:
        s = h.filter((pl.col("cdr3") == c) & (pl.col("v") == V_GENE))
        nas = s.filter(pl.col("state") == "as")["id"].n_unique()
        nh = s.filter(pl.col("state") == "h")["id"].n_unique()
        exact += (nas, nh) == (pas, ph)
        print(f"  {c:<17}{nas:>5}{nh:>4}   {pas:>3}/{ph:<3} {'MATCH' if (nas,nh)==(pas,ph) else ''}",
              flush=True)
    car = h.filter(pl.col("v") == V_GENE).group_by("id").agg(pl.first("state"))
    n_as = car.filter(pl.col("state") == "as").height
    n_h = car.filter(pl.col("state") == "h").height
    print(f"\n  rows matching the paper cell-for-cell: {exact}/8", flush=True)
    print(f"  AS donors carrying >=1 of the 8: {n_as}/25   <-- paper: 15/25 (60%)", flush=True)
    print(f"  healthy donors carrying >=1    : {n_h}/15", flush=True)
    print("  NB the per-clonotype counts run ~1 low on 5 rows: old/ is one BLOOD file per donor,\n"
          "     while the paper's counts pooled blood + sorted + SF samples. The donor union is\n"
          "     what reproduces exactly.", flush=True)


def arm_confirm(clones, donors) -> pl.DataFrame:
    """A: the pre-registered hypothesis. ONE test -> no multiple-testing burden, and well powered
    at 26 vs 12 (where a genome-wide screen would be an honest negative -- cf. covid19 502+/34-)."""
    print("\n=== A: confirmatory -- donor carriage of the 9 (pre-registered) ===", flush=True)
    out = []
    for label, ctrl in (("AS/B27+ vs HD/B27+ (B27-matched)", ("hd", "pos", False)),
                        ("AS/B27+ vs HD/B27+ batch-matched", ("hd", "pos", True)),
                        ("AS/B27+ vs all HD", ("hd", None, False))):
        pos = donors.filter((pl.col("disease_status") == "as") & (pl.col("b27") == "pos"))
        neg = donors.filter(pl.col("disease_status") == ctrl[0])
        if ctrl[1]:
            neg = neg.filter(pl.col("b27") == ctrl[1])
        if ctrl[2]:
            neg = neg.filter(pl.col("proj").is_null())   # drop foreign batches
        for pin in (True, False):
            car = carriage(clones, MOTIF, pin_v=pin)
            a, b, orr, p = _fisher(car, pos["donor"].to_list(), neg["donor"].to_list())
            out.append({"contrast": label, "v_pinned": pin, "n_pos": pos.height, "n_neg": neg.height,
                        "carriers_pos": a, "carriers_neg": b, "odds_ratio": orr, "p_value": p})
            print(f"  {label:36s} V-pinned={str(pin):5s}  {a}/{pos.height} vs {b}/{neg.height}"
                  f"   OR={orr:6.2f}  p={p:.4g}", flush=True)
    print("\n  B27 carriage among HEALTHY donors (expected null if it is disease, not carriage):",
          flush=True)
    hp = donors.filter((pl.col("disease_status") == "hd") & (pl.col("b27") == "pos"))["donor"].to_list()
    hn = donors.filter((pl.col("disease_status") == "hd") & (pl.col("b27") == "neg"))["donor"].to_list()
    a, b, orr, p = _fisher(carriage(clones, MOTIF), hp, hn)
    print(f"    HD/B27+ {a}/{len(hp)} vs HD/B27- {b}/{len(hn)}   OR={orr:.2f}  p={p:.2f}", flush=True)
    print("\n  the AS/B27- donor `Kal` (the internal negative control):", flush=True)
    kal = clones.filter(pl.col("donor") == "Kal")
    print(f"    clones={kal.height:,}  TRBV9/J2-3={kal.filter((pl.col('v_call')==V_GENE) & (pl.col('j_call')==J_GENE)).height}"
          f"  motif hits={kal.filter(pl.col('junction_aa').is_in(MOTIF)).height}", flush=True)
    return pl.DataFrame(out)


def arm_pgen(clones, donors, *, ball: bool = False) -> pl.DataFrame:
    """B: Komech's screen with an exact DP in place of 2e9 Monte-Carlo draws."""
    from vdjtools.model import load_bundled, native
    print("\n=== B: Komech's probabilistic screen, exact DP (they used 2e9 MC draws) ===", flush=True)
    m = load_bundled("TRB", "olga")
    v9 = (clones.filter((pl.col("v_call") == V_GENE) & (pl.col("j_call") == J_GENE)
                        & ~pl.col("junction_aa").str.contains(r"[*~_]"))
          .select("donor", "junction_aa").unique())
    seqs = sorted(v9["junction_aa"].unique().to_list())
    print(f"  unique {V_GENE}/{J_GENE} CDR3s in the cohort: {len(seqs):,}", flush=True)
    with Timer(f"exact aa Pgen x{len(seqs):,}"):
        pg = native.pgen_aa_batch(m, seqs, [V_ALLELE]*len(seqs), [J_ALLELE]*len(seqs), threads=0)
    sharing = v9.group_by("junction_aa").agg(pl.col("donor").n_unique().alias("n_donors"))
    tab = (pl.DataFrame({"junction_aa": seqs, "pgen": pg})
           .join(sharing, on="junction_aa", how="left")
           .with_columns(pl.col("pgen").log10().alias("log10_pgen"),
                         pl.col("junction_aa").is_in(MOTIF).alias("is_motif")))
    if ball:  # 17x the cost of exact and it shifts every member alike -- opt-in, not default
        with Timer(f"1-mismatch ball x{len(seqs):,}"):
            pg1 = native.pgen_aa_batch(m, seqs, [V_ALLELE]*len(seqs), [J_ALLELE]*len(seqs),
                                       mismatches=1, threads=0)
        tab = tab.with_columns(pl.Series("pgen_1mm", pg1))
    # Komech's own statistic (Methods p.1099): per sharing bin, fit a normal to the density of
    # log10 Pgen and flag P<0.001. Their bins are 4, 5, 6, 7 donors with everything ">7 pooled
    # (independent of disease status) to yield >=30 observations per bin".
    #
    # The normal fit was a workaround for MC noise -- with an exact DP the empirical left-tail
    # rank IS the null and needs no distributional assumption. We report BOTH: `p_komech` to show
    # we recover their hits on their terms, `p_empirical` because it is the correct statistic.
    from scipy.stats import norm
    binned = (tab.filter(pl.col("n_donors") >= 4)
              .with_columns(pl.when(pl.col("n_donors") > 7).then(8)
                            .otherwise(pl.col("n_donors")).alias("bin")))  # 8 == ">7, pooled"
    res = []
    for (b,), grp in binned.group_by("bin"):
        x = grp["log10_pgen"].to_numpy()
        if len(x) < 2:   # a bin of one has no dispersion to fit -- report it, do not fake a p
            print(f"  bin {b}: only {len(x)} clonotype(s) -- not scored", flush=True)
            continue
        mu, sd = x.mean(), x.std(ddof=1)
        # LEFT tail: "rarer to generate than its bin-mates, yet shared by this many donors".
        p_norm = norm.cdf(x, mu, sd) if sd > 0 else np.ones_like(x)
        # the same tail, empirically: P(X <= x_i) among the bin's own members. M[i,j] = x_j <= x_i,
        # summed over j.
        p_emp = (x[None, :] <= x[:, None]).sum(1) / len(x)
        res.append(grp.with_columns(pl.Series("p_komech", p_norm),
                                    pl.Series("p_empirical", p_emp)))
    if not res:
        print("  no clonotype shared by >=4 donors", flush=True)
        return tab
    sc = pl.concat(res)
    print(f"\n  sharing bins (paper pools >7): "
          f"{dict(sorted(sc.group_by('bin').len().rows()))}", flush=True)
    sig = sc.filter(pl.col("p_komech") < 0.001)
    n_test = sc.filter(pl.col("is_motif")).height
    print(f"  shared by >=4 donors: {sc.height}   below Komech's P<0.001: {sig.height}", flush=True)
    print(f"  of those, motif members: {sig.filter(pl.col('is_motif')).height} of {n_test} testable"
          f"  (enrichment: {sig.filter(pl.col('is_motif')).height}/{sig.height} of the hits are "
          f"motif, vs {n_test}/{sc.height} of the tested)", flush=True)
    print(f"\n  {'CDR3':<17}{'donors':>7}{'log10 Pgen':>12}{'p_komech':>10}{'p_emp':>8}{'motif':>7}",
          flush=True)
    for r in sc.sort("p_komech").head(12).iter_rows(named=True):
        print(f"  {r['junction_aa']:<17}{r['n_donors']:>7}{r['log10_pgen']:>12.2f}"
              f"{r['p_komech']:>10.2e}{r['p_empirical']:>8.3f}"
              f"{'  YES' if r['is_motif'] else '':>7}", flush=True)
    mm = tab.filter(pl.col("is_motif"))
    if ball and mm.height:
        shift = np.log10(mm["pgen_1mm"].to_numpy() / mm["pgen"].to_numpy())
        print(f"\n  1mm ball vs exact, over the motif: median {np.median(shift):.2f} log10, "
              f"spread {shift.max()-shift.min():.2f}", flush=True)
        print("  => it lifts every member alike; a covariate, not a discriminator.", flush=True)
    return sc


def arm_screen(clones, donors, *, min_incidence: int = 4) -> pl.DataFrame:
    """C: discovery -- find it WITHOUT being told the sequences.

    Two levers make BH viable at 26 vs 12 (where covid19's 502+/34- gave 0 on 52k features):
    V/J-pinning cuts the feature space ~220x, and an incidence floor drops the singleton tail that
    can never reach significance anyway. `min_incidence=4` is not a tuned knob -- it is Komech's
    own screening criterion ("shared between at least any four donors", p.1099).
    """
    print("\n=== C: V/J-pinned discovery (no knowledge of the 9) ===", flush=True)
    pos = set(donors.filter((pl.col("disease_status") == "as") & (pl.col("b27") == "pos"))["donor"])
    neg = set(donors.filter((pl.col("disease_status") == "hd") & (pl.col("b27") == "pos"))["donor"])
    keep = pos | neg
    sub = clones.filter(pl.col("donor").is_in(list(keep))
                        & ~pl.col("junction_aa").str.contains(r"[*~_]"))
    for label, f in (("unpinned (junction_aa only)", pl.lit(True)),
                     (f"pinned to {V_GENE}/{J_GENE}", (pl.col("v_call") == V_GENE)
                      & (pl.col("j_call") == J_GENE))):
        d = sub.filter(f).select("donor", "junction_aa").unique()
        n_feat = d["junction_aa"].n_unique()
        print(f"  {label:32s} features={n_feat:>7,}   min attainable BH p "
              f"~ {0.05/n_feat:.1e}", flush=True)
    d = (sub.filter((pl.col("v_call") == V_GENE) & (pl.col("j_call") == J_GENE))
         .select("donor", "junction_aa").unique())
    inc = d.group_by("junction_aa").agg(
        pl.col("donor").is_in(list(pos)).sum().alias("a"),
        pl.col("donor").is_in(list(neg)).sum().alias("b"))
    n_all = inc.height
    inc = inc.filter((pl.col("a") + pl.col("b")) >= min_incidence)
    print(f"  incidence floor >={min_incidence} donors (Komech's own criterion): "
          f"{n_all:,} -> {inc.height:,} features", flush=True)
    a = inc["a"].to_numpy(); b = inc["b"].to_numpy()
    c = len(pos) - a; dd = len(neg) - b
    p = bstats.fisher_p(a, b, c, dd, alternative="greater")
    # Ties matter here: Fisher on small integer counts produces long ties, and polars' group_by
    # order is not deterministic -- sorting on p_value alone made the motif's ranks jump between
    # identical runs (1,2,4,7,10,35,43,122 vs 1,2,5,7,10,21,22,120). Break ties on the key so the
    # ranking is reproducible. (Same failure SOURCES.md records for `select_candidates`.)
    out = (inc.with_columns(pl.Series("p_value", p),
                            pl.Series("q_value", bstats.fdr_bh(p)),
                            pl.Series("odds_ratio", bstats.odds_ratio(a, b, c, dd)))
           .with_columns(pl.col("junction_aa").is_in(MOTIF).alias("is_motif"))
           .sort(["p_value", "junction_aa"]))
    sig = out.filter(pl.col("q_value") < 0.05)
    print(f"\n  tested {out.height:,} pinned features; BH q<0.05: {sig.height}", flush=True)
    print(f"  {'rank':>5}  {'CDR3':<17}{'AS':>4}{'HD':>4}{'OR':>8}{'p':>10}{'q':>9}{'motif':>7}",
          flush=True)
    for i, r in enumerate(out.head(10).iter_rows(named=True), 1):
        print(f"  {i:>5}  {r['junction_aa']:<17}{r['a']:>4}{r['b']:>4}{r['odds_ratio']:>8.1f}"
              f"{r['p_value']:>10.2e}{r['q_value']:>9.3f}{'  YES' if r['is_motif'] else '':>7}",
              flush=True)
    ranks = [i for i, r in enumerate(out.iter_rows(named=True), 1) if r["is_motif"]]
    print(f"\n  motif members' ranks among {out.height:,}: {ranks}", flush=True)

    # BH over 273 features at 26-vs-12 cannot clear 0.05 -- the smallest attainable p is Fisher's
    # 12/26-vs-0/12 = 3.6e-3, and 3.6e-3 * 273 > 0.05. That is the covid19 lesson (502+/34- gave 0
    # on 52k features), not a property of the motif. The powered question is not "does any single
    # clonotype survive BH" but "is the FAMILY concentrated at the top of the ranking" -- one
    # hypergeometric over the ranking, no per-feature burden.
    from scipy.stats import hypergeom
    n_motif = len(ranks)
    print(f"\n  is the family concentrated at the top? (hypergeometric over the ranking)", flush=True)
    print(f"  {'top-k':>6}{'motif in top-k':>16}{'expected':>10}{'p':>12}", flush=True)
    for k in (10, 25, 50):
        obs = sum(r <= k for r in ranks)
        exp = k * n_motif / out.height
        p = hypergeom.sf(obs - 1, out.height, n_motif, k)
        print(f"  {k:>6}{obs:>16}{exp:>10.2f}{p:>12.2e}", flush=True)
    return out


def arm_oracle(screen: pl.DataFrame, vdjdb: Path) -> None:
    """VDJdb cross-check. NB this oracle is PARTLY CIRCULAR: Yang 2022 (Nature 612:771) is the
    sequel to Komech 2018 by the same group, and its TCRs came from AS patients -- plausibly these
    donors. It confirms restriction+epitope; it does NOT independently validate discovery."""
    print("\n=== oracle: VDJdb (Yang 2022 Nature -- same group; partly circular) ===", flush=True)
    if not vdjdb.exists():
        print(f"  SKIP -- no VDJdb at {vdjdb}", flush=True)
        return
    db = pl.read_csv(vdjdb, separator="\t", infer_schema_length=0)
    hit = db.filter(pl.col("cdr3").is_in(MOTIF + [MOTIF_VDJDB]))
    print(f"  exact VDJdb records for the 9 (+the 10th): {hit.height}", flush=True)
    if hit.height:
        print(hit.select("cdr3", "v.segm", "j.segm", "antigen.epitope", "antigen.gene",
                         "antigen.species", "mhc.a", "reference.id").unique()
              .sort("cdr3"), flush=True)
    # the germline trap: the J2-3 ending alone is NOT the motif.
    fam = db.filter(pl.col("cdr3").str.ends_with("STDTQYF") & pl.col("v.segm").str.contains(V_GENE))
    print(f"\n  {V_GENE} + *STDTQYF in VDJdb: {fam.height} records, restricted by "
          f"{sorted(set(fam['mhc.a'].to_list()))[:4]}...", flush=True)
    print("  => STDTQYF is germline TRBJ2-3, not the motif. Matching the ending is an FP "
          "generator; the specificity is in the central residues.", flush=True)
    if MOTIF_VDJDB in screen["junction_aa"].to_list():
        r = screen.filter(pl.col("junction_aa") == MOTIF_VDJDB).row(0, named=True)
        print(f"\n  PROSPECTIVE: {MOTIF_VDJDB} (in VDJdb, NOT in Komech's 9) -> "
              f"AS {r['a']} / HD {r['b']}, OR={r['odds_ratio']:.1f}, p={r['p_value']:.3g}",
              flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--arm", default="all",
                    choices=["all", "cohort", "fig1", "confirm", "pgen", "screen"])
    ap.add_argument("--vdjdb", type=Path, default=Path(DEF_VDJDB))
    ap.add_argument("--ball", action="store_true",
                    help="also compute the 1-mismatch Pgen ball (~17x slower)")
    ap.add_argument("--min-incidence", type=int, default=4,
                    help="donor floor for the discovery screen (Komech used 4)")
    ap.add_argument("--out", type=Path, default=RES)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    snap = _snapshot()
    with Timer("load new/ (68 samples)"):
        clones, donors = load_new(snap)
    print(f"  {clones.height:,} rows, {donors.height} donors", flush=True)

    if args.arm in ("all", "cohort"):
        arm_cohort(clones, donors)
    if args.arm in ("all", "fig1"):
        with Timer("old/ Fig.1 reproduction"):
            arm_fig1(snap)
    if args.arm in ("all", "confirm"):
        arm_confirm(clones, donors).write_parquet(args.out / "as_b27_confirm.parquet")
    if args.arm in ("all", "pgen"):
        arm_pgen(clones, donors, ball=args.ball).write_parquet(
            args.out / "as_b27_pgen.parquet")
    if args.arm in ("all", "screen"):
        sc = arm_screen(clones, donors, min_incidence=args.min_incidence)
        sc.write_parquet(args.out / "as_b27_screen.parquet")
        arm_oracle(sc, args.vdjdb)

    print(f"\nTotal {time.perf_counter()-t0:.0f}s   peak RSS {_rss_gb():.2f} GB", flush=True)
    print(f"-> {args.out}", flush=True)


if __name__ == "__main__":
    main()
