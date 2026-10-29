"""Independent productive-input software validation of isalgo and frozen OLGA models.

Fresh functional reads are outside the all-nonfunctional fitting inputs. They are
selection-shaped observations, not unbiased recombination heldout likelihood.
The accepted protocol remains first20000 raw reads, then first N_CMP=300 eligible
unique junctions. No annotation/parquet cache or unproven OOF holdout is used.

BUILD_GROUPS=human,mouse selects species; LOCI optionally restricts loci.
RAW_FASTQ_DIR points to the existing raw mirror (<group>/<locus>.functional.fq.gz).
VDJTOOLS_BUNDLE_ROOT selects the bundle parent containing isalgo/<context>.
VDJTOOLS_OLGA_MODELS overrides human oracle fixtures; VDJTOOLS_MOUSE_OLGA_MODELS
optionally overrides OLGA-pip's mouse default_models. EM_THREADS bounds native work.
Without RAW_FASTQ_DIR, the existing authorized data.fetch_fastq route is used.
EM_WORK optionally selects the parent of fresh temporary annotation directories.
"""
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import olga
import olga.generation_probability as ogp
import olga.load_model as olm
import polars as pl

from vdjtools.model import Model, data, native
from vdjtools.model.io import from_olga
from vdjtools.model.schema import normalization_keys

BUNDLED = Path(os.environ.get("VDJTOOLS_BUNDLE_ROOT", "python/vdjtools/model/_bundled"))
# TRG/TRD fixtures originate from mirpy legacy-v2 commit aeccd75 and were verified
# byte-identical to the files used to build the frozen bundled OLGA parquet.
_REPO_OLGA = Path(__file__).resolve().parent.parent / "tests/python/fixtures/olga/default_models"
OLGA = Path(os.environ.get("VDJTOOLS_OLGA_MODELS", str(_REPO_OLGA)))
MOUSE_OLGA = Path(os.environ.get("VDJTOOLS_MOUSE_OLGA_MODELS", str(Path(olga.__file__).parent / "default_models")))
WORK = Path(os.environ.get("EM_WORK", "/tmp/em_work"))
N = int(os.environ.get("N_CMP", "300"))
THREADS = int(os.environ.get("EM_THREADS", "0"))
LOCI = {"TRA": ("T_alpha", "VJ"), "TRB": ("T_beta", "VDJ"),
        "TRG": ("T_gamma", "VJ"), "TRD": ("T_delta", "VDJ"),
        "IGH": ("B_heavy", "VDJ"), "IGK": ("B_kappa", "VJ"), "IGL": ("B_lambda", "VJ")}
CONTEXTS = {"human": tuple(LOCI), "mouse": ("TRA", "TRB", "IGH", "IGK", "IGL")}


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


def _gene_conditioned_pgens(model, seqs, v_calls, j_calls, *, threads=THREADS):
    """Sum disjoint explicit gene/allele states; unknown-only calls have zero support."""
    from vdjtools.model.infer import _gene_to_alleles, call_alleles

    if not len(seqs) == len(v_calls) == len(j_calls):
        raise ValueError("sequences and V/J calls must have equal lengths")
    vi, ji = _gene_to_alleles(model, "v"), _gene_to_alleles(model, "j")
    expanded, vs, js, owners = [], [], [], []
    for i, (sequence, v, j) in enumerate(zip(seqs, v_calls, j_calls)):
        for vv in call_alleles(vi, v):
            for jj in call_alleles(ji, j):
                expanded.append(sequence.upper())
                vs.append(vv)
                js.append(jj)
                owners.append(i)
    result = np.zeros(len(seqs), dtype=float)
    if expanded:
        scores = np.asarray(native.pgen_nt_batch(model, expanded, vs, js, threads=threads))
        np.add.at(result, owners, scores)
    return result


def held_out_loglik(model, seqs, v_calls, j_calls):
    """Compatibility helper: positive-only log10 mean, zero count and total count.

    This function does not establish provenance or independence of supplied observations.
    """
    summary = _pgen_summary(_gene_conditioned_pgens(model, seqs, v_calls, j_calls))
    return summary["mean_log10_positive"], summary["zero"], summary["total"]


def _pgen_summary(values):
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or (values > 1 + 1e-9).any():
        raise ValueError("invalid Pgen values")
    positive = values > 0
    total, nonzero = len(values), int(positive.sum())
    positive_mean = float(np.log10(values[positive]).mean()) if nonzero else None
    return dict(total=total, positive=nonzero, zero=total - nonzero,
                mean_pgen=float(values.mean()) if total else None,
                mean_log10_positive=positive_mean,
                mean_log10_all=positive_mean if nonzero == total else None,
                mean_log10_all_is_negative_infinity=bool(total and nonzero < total))


def functional_queries(group, locus):
    """Fresh accepted functional-input protocol; no read identifiers leave this helper."""
    if N < 1:
        raise ValueError("N_CMP must be positive")
    raw_dir = os.environ.get("RAW_FASTQ_DIR")
    raw = Path(raw_dir) / group / f"{locus}.functional.fq.gz" if raw_dir else data.fetch_fastq(group, locus, "functional")
    WORK.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=f"validate-{group}-{locus}-", dir=WORK) as scratch:
        reads = data.annotate_reads(raw, out_dir=scratch, prefix=f"{group}_{locus}_functional",
                                    organism=group, cap=20000, threads=THREADS)
        uniq = data.unique_clonotypes(reads).filter(
            pl.col("junction").str.to_uppercase().str.contains(r"^[ACGT]+$")
            & ((pl.col("junction").str.len_bytes() % 3) == 0))
        test = uniq.head(N)
        if not test.height:
            raise ValueError(f"no eligible functional junctions for {group}/{locus}")
        counts = dict(raw_read_cap=20000, annotated_reads=reads.height,
                      eligible_clonotypes=uniq.height, requested_clonotypes=N, selected_clonotypes=test.height)
        return (test["junction"].str.to_uppercase().to_list(), test["v_call"].to_list(),
                test["j_call"].to_list(), counts)


def run(locus: str, group: str = "human") -> dict:
    name, chain = LOCI[locus]
    oracle_path = (OLGA if group == "human" else MOUSE_OLGA) / f"{group}_{name}"
    key = locus if group == "human" else f"mouse_{locus}"
    # Each model retains its own reference; ARDA and OLGA germlines are never mixed.
    m_olga = (Model.load(BUNDLED / "olga" / locus) if group == "human" else
              from_olga(oracle_path, locus=locus, organism=group))
    baseline = "olga-legacy" if group == "human" else "olga-pip-import"
    m_isalgo = Model.load(BUNDLED / "isalgo" / key)
    opg = olga_pgen_model(oracle_path, chain)
    seqs, v_calls, j_calls, counts = functional_queries(group, locus)
    print(f"\n{group}/{locus}: independent selected productive inputs; {counts}", flush=True)
    ours = np.asarray(native.pgen_nt_batch(m_olga, seqs, threads=THREADS))
    theirs = np.asarray([opg.compute_nt_CDR3_pgen(sequence) for sequence in seqs])
    # Include every zero in the oracle assertion; filtering both-positive vectors hides defects.
    np.testing.assert_allclose(ours, theirs, rtol=1e-9, atol=0,
                               err_msg=f"{group}/{locus}: exact OLGA oracle mismatch")
    both_positive = (ours > 0) & (theirs > 0)
    positive_ours, positive_theirs = np.log10(ours[both_positive]), np.log10(theirs[both_positive])
    correlation = (float(np.corrcoef(positive_ours, positive_theirs)[0, 1])
                   if both_positive.sum() >= 2 and np.ptp(positive_ours) > 0 and np.ptp(positive_theirs) > 0
                   else None)
    max_relative_error = (float(np.max(np.abs(ours[both_positive] - theirs[both_positive]) /
                                       theirs[both_positive])) if both_positive.any() else None)
    oracle_summary = dict(total=len(seqs), both_positive=int(both_positive.sum()),
                          max_relative_error_both_positive=max_relative_error,
                          log10_correlation_both_positive=correlation, ours_zero=int((ours == 0).sum()),
                          oracle_zero=int((theirs == 0).sum()), rtol=1e-9, atol=0,
                          all_values_agree=True)
    print(f"ORACLE full-vector agreement: {oracle_summary}")
    tv_rows = []
    print("COMMON EVENTS: descriptive per-parent TV; only parents with mass in both models")
    for ev in ("v_3_del", "j_5_del", "d_del", "vd_ins", "dj_ins", "vd_dinucl", "dj_dinucl",
               "d_gene", "j_choice", "vj_ins", "vj_dinucl"):
        if ev not in m_olga.tables or ev not in m_isalgo.tables:
            continue
        a, b = m_olga.tables[ev], m_isalgo.tables[ev]
        parents = list(normalization_keys(m_olga.manifest.events[ev]))
        tvs = _per_parent_tv(a, b, parents)
        if not tvs:
            continue
        row = dict(event=ev, mean_tv=float(np.mean(tvs)), n_parents=len(tvs),
                   entropy_olga=_entropy(a["p"].to_numpy()), entropy_isalgo=_entropy(b["p"].to_numpy()))
        tv_rows.append(row)
        print(row)
    usage = {label: v_support(model) for label, model in ((baseline, m_olga), ("isalgo", m_isalgo))}
    print(f"V GENE USAGE: descriptive, protocol-dependent: {usage}")
    pgens = {}
    for label, model in ((baseline, m_olga), ("isalgo", m_isalgo)):
        pgens[label] = _pgen_summary(_gene_conditioned_pgens(model, seqs, v_calls, j_calls))
        print(f"PRODUCTIVE-INPUT GENE-CALL PGEN {label}: {pgens[label]}")
    insertion_zero = {
        label: {event: float(table.filter(pl.col("length") == 0)["p"].sum())
                for event, table in model.tables.items()
                if event in ("vj_ins", "vd_ins", "dj_ins", "dd_ins")}
        for label, model in ((baseline, m_olga), ("isalgo", m_isalgo))}
    print(f"INSERTION LENGTH P(N=0), conditional on using each event: {insertion_zero}")
    return dict(group=group, locus=locus, population="independent selected productive inputs",
                baseline=baseline, queries=counts, oracle=oracle_summary, common_event_tv=tv_rows,
                v_usage=usage, pgen=pgens, insertion_zero=insertion_zero)


def main():
    groups = os.environ.get("BUILD_GROUPS", "human,mouse").split(",")
    if any(group not in CONTEXTS for group in groups):
        raise ValueError("BUILD_GROUPS must select human and/or mouse")
    want = os.environ.get("LOCI")
    if want and set(want.split(",")) - set(LOCI):
        raise ValueError("LOCI contains an unknown locus")
    reports = []
    for group in groups:
        loci = tuple(want.split(",")) if want else CONTEXTS[group]
        if set(loci) - set(CONTEXTS[group]):
            raise ValueError(f"requested loci are unsupported for {group}")
        reports.extend(run(locus, group) for locus in loci)
    print(json.dumps(dict(contexts=reports), allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
