"""Model I/O: import from OLGA's format, and save/load our native parquet directory.

``from_olga`` uses OLGA's own parser (an optional, import-time dependency — the ``[oracle]``
extra) to read a bootstrap model, then tabularises its processed arrays into our long-format
polars schema. Native ``load_model`` reads the parquet directory and needs no OLGA.

Conventions reproduced exactly (so a loaded model is bit-faithful to OLGA — verified by the
round-trip test):

- Gene order follows OLGA's IGoR gene index (``genV``/``genD``/``genJ`` list order).
- ``ndel = array_index - max_palindrome`` (biological deletion; negatives are palindromic P-nt).
- ``PDJ`` (VDJ) is factored to ``j_choice`` (marginal ``P(J)``) × ``d_gene`` given J (``P(D|J)``);
  ``PVJ`` (VJ) to ``v_choice`` (``P(V)``) × ``j_choice`` given V (``P(J|V)``). Both reconstruct exactly.
- Dinucleotide row ``(from_nt, to_nt, p)`` = ``R[to_nt, from_nt]`` (OLGA's column-stochastic Markov).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

from .events import Event, EventKind
from .model import Model
from .schema import Manifest

_NT = "ACGT"


def _gene(allele: str) -> str:
    return allele.split("*")[0]


def _read_anchors(path: Path) -> dict[str, tuple[int, str]]:
    """Anchor CSV -> {gene_allele: (anchor_index, functionality_char)}."""
    if not path.exists():
        return {}
    df = pl.read_csv(path)
    return {r[0]: (int(r[1]), str(r[2])) for r in df.iter_rows()}


def _genomic_table(gen_list, anchors, cut_segs, *, seg: str, has_anchor: bool) -> pl.DataFrame:
    """Build a genes_<seg> frame from OLGA's gen* list + cut segments (+ anchors for V/J)."""
    rows = []
    for i, entry in enumerate(gen_list):
        name = entry[0]
        if has_anchor:  # genV/genJ entries are [name, cdr3_trim, full_germline]
            cdr3_seg, full = entry[1], entry[2]
            anchor = anchors.get(name, (-1, ""))[0]
        else:  # genD entries are [name, germline]
            full, cdr3_seg, anchor = entry[1], entry[1], -1
        rows.append(
            {
                f"{seg}_allele": name,
                "gene": _gene(name),
                "full_germline": full,
                "cdr3_segment": cdr3_seg,
                "cut_segment": cut_segs[i],
                "anchor": anchor,
                "functional": len(cut_segs[i]) > 0,
            }
        )
    return pl.DataFrame(rows)


def _gene_choice_table(alleles: list[str], p: np.ndarray, col: str) -> pl.DataFrame:
    return pl.DataFrame({col: alleles, "p": p.astype(np.float64)})


def _deletion_table(alleles: list[str], p_del_given: np.ndarray, max_pal: int, col: str) -> pl.DataFrame:
    """p_del_given is [del_index, gene]; emit (gene_allele, ndel, p) with ndel = index - max_pal."""
    n_idx, n_gene = p_del_given.shape
    ndel = np.arange(n_idx) - max_pal
    a = np.repeat(alleles, n_idx)
    d = np.tile(ndel, n_gene)
    return pl.DataFrame(
        {col: a, "ndel": d.astype(np.int16), "p": p_del_given.T.reshape(-1).astype(np.float64)}
    )


def _ins_len_table(p_ins: np.ndarray) -> pl.DataFrame:
    return pl.DataFrame(
        {"length": np.arange(len(p_ins), dtype=np.int16), "p": p_ins.astype(np.float64)}
    )


def _dinucl_table(R: np.ndarray) -> pl.DataFrame:
    """R[next, prev], columns sum to 1 -> rows (from_nt=prev, to_nt=next, p=P(next|prev))."""
    frm = np.repeat(np.arange(4), 4)  # prev
    to = np.tile(np.arange(4), 4)  # next
    p = R[to, frm]
    return pl.DataFrame(
        {"from_nt": frm.astype(np.uint8), "to_nt": to.astype(np.uint8), "p": p.astype(np.float64)}
    )


def from_olga(model_dir: str | Path, *, locus: str, organism: str = "human") -> Model:
    """Import an OLGA default-model directory into a :class:`Model`.

    Args:
        model_dir: Directory with ``model_params.txt``, ``model_marginals.txt`` and the
            ``V/J_gene_CDR3_anchors.csv`` files (e.g. an OLGA ``default_models/*`` folder).
        locus: Locus label to record (e.g. ``"TRB"``).
        organism: Organism label (default ``"human"``).

    Returns:
        A validated :class:`Model` whose tables reproduce OLGA's arrays exactly.
    """
    import olga.load_model as olm  # optional dep (the [oracle] extra)

    d = Path(model_dir)
    params, marg = d / "model_params.txt", d / "model_marginals.txt"
    v_anchors = _read_anchors(d / "V_gene_CDR3_anchors.csv")
    j_anchors = _read_anchors(d / "J_gene_CDR3_anchors.csv")
    chain_type = "VDJ" if "@d_gene" in marg.read_text() else "VJ"

    tables: dict[str, pl.DataFrame] = {}
    events: dict[str, Event] = {}

    if chain_type == "VDJ":
        g = olm.GenomicDataVDJ()
        g.load_igor_genomic_data(str(params), str(d / "V_gene_CDR3_anchors.csv"), str(d / "J_gene_CDR3_anchors.csv"))
        m = olm.GenerativeModelVDJ()
        m.load_and_process_igor_model(str(marg))

        v_alleles = [x[0] for x in g.genV]
        d_alleles = [x[0] for x in g.genD]
        j_alleles = [x[0] for x in g.genJ]

        pj = m.PDJ.sum(axis=0)  # P(J), shape (J,)
        # P(D | J): PDJ[d, j] / P(J); where P(J)=0, use a uniform placeholder (reconstruction is 0 either way).
        safe_pj = np.where(pj > 0, pj, 1.0)
        p_d_given_j = m.PDJ / safe_pj[np.newaxis, :]
        p_d_given_j[:, pj == 0] = 1.0 / len(d_alleles)

        tables["v_choice"] = _gene_choice_table(v_alleles, m.PV, "v_allele")
        tables["j_choice"] = _gene_choice_table(j_alleles, pj, "j_allele")
        tables["d_gene"] = pl.DataFrame(
            {
                "j_allele": np.repeat(j_alleles, len(d_alleles)),
                "d_allele": np.tile(d_alleles, len(j_alleles)),
                "p": p_d_given_j.T.reshape(-1).astype(np.float64),  # [j, d] order
            }
        )
        tables["n_d"] = pl.DataFrame({"n_d": np.array([1], np.uint8), "p": np.array([1.0])})
        tables["v_3_del"] = _deletion_table(v_alleles, m.PdelV_given_V, g.max_delV_palindrome, "v_allele")
        tables["j_5_del"] = _deletion_table(j_alleles, m.PdelJ_given_J, g.max_delJ_palindrome, "j_allele")
        tables["d_del"] = _d_del_table(d_alleles, m.PdelDldelDr_given_D, g.max_delDl_palindrome, g.max_delDr_palindrome)
        tables["vd_ins"] = _ins_len_table(m.PinsVD)
        tables["dj_ins"] = _ins_len_table(m.PinsDJ)
        tables["vd_dinucl"] = _dinucl_table(m.Rvd)
        tables["dj_dinucl"] = _dinucl_table(m.Rdj)

        events = {
            "v_choice": Event("v_choice", EventKind.GENE_CHOICE),
            "j_choice": Event("j_choice", EventKind.GENE_CHOICE),
            "d_gene": Event("d_gene", EventKind.GENE_CHOICE, ("j_choice",)),
            "n_d": Event("n_d", EventKind.N_D),
            "v_3_del": Event("v_3_del", EventKind.DELETION, ("v_choice",)),
            "j_5_del": Event("j_5_del", EventKind.DELETION, ("j_choice",)),
            "d_del": Event("d_del", EventKind.DELETION_2D, ("d_gene",)),
            "vd_ins": Event("vd_ins", EventKind.INS_LENGTH),
            "dj_ins": Event("dj_ins", EventKind.INS_LENGTH),
            "vd_dinucl": Event("vd_dinucl", EventKind.DINUCLEOTIDE),
            "dj_dinucl": Event("dj_dinucl", EventKind.DINUCLEOTIDE),
        }
        palindrome_max = {
            "v_3": int(g.max_delV_palindrome),
            "d_5": int(g.max_delDl_palindrome),
            "d_3": int(g.max_delDr_palindrome),
            "j_5": int(g.max_delJ_palindrome),
        }
        genomic = {
            "genes_v": _genomic_table(g.genV, v_anchors, g.cutV_genomic_CDR3_segs, seg="v", has_anchor=True),
            "genes_d": _genomic_table(g.genD, {}, g.cutD_genomic_CDR3_segs, seg="d", has_anchor=False),
            "genes_j": _genomic_table(g.genJ, j_anchors, g.cutJ_genomic_CDR3_segs, seg="j", has_anchor=True),
        }
    else:  # VJ
        g = olm.GenomicDataVJ()
        g.load_igor_genomic_data(str(params), str(d / "V_gene_CDR3_anchors.csv"), str(d / "J_gene_CDR3_anchors.csv"))
        m = olm.GenerativeModelVJ()
        m.load_and_process_igor_model(str(marg))

        v_alleles = [x[0] for x in g.genV]
        j_alleles = [x[0] for x in g.genJ]

        pv = m.PVJ.sum(axis=1)  # P(V), shape (V,)
        safe_pv = np.where(pv > 0, pv, 1.0)
        p_j_given_v = m.PVJ / safe_pv[:, np.newaxis]
        p_j_given_v[pv == 0, :] = 1.0 / len(j_alleles)

        tables["v_choice"] = _gene_choice_table(v_alleles, pv, "v_allele")
        tables["j_choice"] = pl.DataFrame(
            {
                "v_allele": np.repeat(v_alleles, len(j_alleles)),
                "j_allele": np.tile(j_alleles, len(v_alleles)),
                "p": p_j_given_v.reshape(-1).astype(np.float64),  # [v, j] order
            }
        )
        tables["v_3_del"] = _deletion_table(v_alleles, m.PdelV_given_V, g.max_delV_palindrome, "v_allele")
        tables["j_5_del"] = _deletion_table(j_alleles, m.PdelJ_given_J, g.max_delJ_palindrome, "j_allele")
        tables["vj_ins"] = _ins_len_table(m.PinsVJ)
        tables["vj_dinucl"] = _dinucl_table(m.Rvj)

        events = {
            "v_choice": Event("v_choice", EventKind.GENE_CHOICE),
            "j_choice": Event("j_choice", EventKind.GENE_CHOICE, ("v_choice",)),
            "v_3_del": Event("v_3_del", EventKind.DELETION, ("v_choice",)),
            "j_5_del": Event("j_5_del", EventKind.DELETION, ("j_choice",)),
            "vj_ins": Event("vj_ins", EventKind.INS_LENGTH),
            "vj_dinucl": Event("vj_dinucl", EventKind.DINUCLEOTIDE),
        }
        palindrome_max = {"v_3": int(g.max_delV_palindrome), "j_5": int(g.max_delJ_palindrome)}
        genomic = {
            "genes_v": _genomic_table(g.genV, v_anchors, g.cutV_genomic_CDR3_segs, seg="v", has_anchor=True),
            "genes_j": _genomic_table(g.genJ, j_anchors, g.cutJ_genomic_CDR3_segs, seg="j", has_anchor=True),
        }

    manifest = Manifest(
        locus=locus,
        organism=organism,
        chain_type=chain_type,
        events=events,
        palindrome_max=palindrome_max,
        source=f"olga:{d.name}",
    )
    return Model(manifest=manifest, tables=tables, genomic=genomic).validate()


def _d_del_table(d_alleles, p_joint: np.ndarray, max_dl: int, max_dr: int) -> pl.DataFrame:
    """PdelDldelDr_given_D is [delDl, delDr, D]; emit (d_allele, ndel5, ndel3, p)."""
    n5, n3, nd = p_joint.shape
    ndel5 = np.arange(n5) - max_dl
    ndel3 = np.arange(n3) - max_dr
    # index order: d (outer), then ndel5, then ndel3 (inner)
    a = np.repeat(d_alleles, n5 * n3)
    g5 = np.tile(np.repeat(ndel5, n3), nd)
    g3 = np.tile(ndel3, n5 * nd)
    p = np.transpose(p_joint, (2, 0, 1)).reshape(-1)  # [D, delDl, delDr] -> flat
    return pl.DataFrame(
        {
            "d_allele": a,
            "ndel5": g5.astype(np.int16),
            "ndel3": g3.astype(np.int16),
            "p": p.astype(np.float64),
        }
    )


# --------------------------------------------------------------------------------------------
# arda-native model construction (germline + names from arda; marginals are placeholders that
# EM (`infer_native`) refits). Generation/Pgen read only `{seg}_allele`, `cut_segment`,
# `functional` from `genomic`, so an arda-built genomic yields arda IMGT allele names verbatim.
# --------------------------------------------------------------------------------------------

_DEFAULT_PALINDROME_MAX = {"v_3": 4, "d_5": 4, "d_3": 4, "j_5": 4}


def _decay_p(values: np.ndarray, scale: float) -> np.ndarray:
    """Normalized weights ``∝ exp(-|value|/scale)`` — a placeholder favouring small trims/inserts
    while keeping full support so EM can move the mass."""
    w = np.exp(-np.abs(values) / scale)
    return w / w.sum()


def _arda_vj_genomic(gl_seg: pl.DataFrame, full: dict, seg: str, max_pal: int) -> pl.DataFrame:
    """genes_v / genes_j frame from arda: cut_segment via palindrome extension; full_germline +
    anchor from arda scaffolds (``""``/``-1`` where absent); functional = F & has cut_segment."""
    from . import reference as ref

    S = seg.upper()
    rows = []
    for r in gl_seg.iter_rows(named=True):
        allele, cdr3 = r["allele"], r["sequence"]
        cut = ref.cut_segment(cdr3, S, max_pal) if cdr3 else ""
        if set(cut) - set("ACGT"):  # ambiguous (N) germline — unusable in the native DP; drop
            continue
        fg, anchor = full.get((S, allele), ("", -1))
        rows.append({
            f"{seg}_allele": allele, "gene": r["gene"], "full_germline": fg,
            "cdr3_segment": cdr3, "cut_segment": cut, "anchor": anchor,
            "functional": bool(r["functional"]) and len(cut) > 0,
        })
    return pl.DataFrame(rows)


def _arda_d_genomic(gl_d: pl.DataFrame, max5: int, max3: int) -> pl.DataFrame:
    """genes_d frame from arda D germline (full = CDR3-region = D; anchor -1)."""
    from . import reference as ref

    rows = []
    for r in gl_d.iter_rows(named=True):
        seq = r["sequence"]
        cut = ref.cut_segment_d(seq, max5, max3)
        if set(cut) - set("ACGT"):  # ambiguous (N) germline — unusable in the native DP; drop
            continue
        rows.append({
            "d_allele": r["allele"], "gene": r["gene"], "full_germline": seq,
            "cdr3_segment": seq, "cut_segment": cut, "anchor": -1, "functional": len(cut) > 0,
        })
    return pl.DataFrame(rows)


def _uniform_choice(alleles: list[str], functional: list[bool], col: str) -> pl.DataFrame:
    """Uniform gene choice over *functional* alleles (0 for the rest); sums to 1."""
    f = np.asarray(functional, dtype=bool)
    n = int(f.sum())
    p = (f / n) if n else np.full(len(alleles), 1.0 / len(alleles))
    return pl.DataFrame({col: alleles, "p": p.astype(np.float64)})


def _uniform_cond_choice(parents: list[str], children: list[str], child_functional: list[bool],
                         parent_col: str, child_col: str) -> pl.DataFrame:
    """P(child|parent) uniform over functional children, identical for every parent."""
    f = np.asarray(child_functional, dtype=bool)
    n = int(f.sum())
    pc = (f / n) if n else np.full(len(children), 1.0 / len(children))
    return pl.DataFrame({
        parent_col: np.repeat(parents, len(children)),
        child_col: np.tile(children, len(parents)),
        "p": np.tile(pc, len(parents)).astype(np.float64),
    })


def _placeholder_deletion(alleles: list[str], ndel_min: int, ndel_max: int, col: str) -> pl.DataFrame:
    """Per-allele deletion table ``∝ exp(-|ndel|/3)`` over ``[ndel_min, ndel_max]`` (wide support)."""
    ndel = np.arange(ndel_min, ndel_max + 1)
    p = _decay_p(ndel.astype(float), 3.0)
    return pl.DataFrame({
        col: np.repeat(alleles, len(ndel)),
        "ndel": np.tile(ndel, len(alleles)).astype(np.int16),
        "p": np.tile(p, len(alleles)).astype(np.float64),
    })


def _placeholder_d_del(d_alleles: list[str], n5min: int, n5max: int, n3min: int, n3max: int) -> pl.DataFrame:
    """Per-D joint 5'/3' deletion ``∝ exp(-(|n5|+|n3|)/3)``, normalized per D."""
    ndel5, ndel3 = np.arange(n5min, n5max + 1), np.arange(n3min, n3max + 1)
    w = np.exp(-np.abs(ndel5)[:, None] / 3.0) * np.exp(-np.abs(ndel3)[None, :] / 3.0)
    p = (w / w.sum()).reshape(-1)
    n5, n3 = len(ndel5), len(ndel3)
    return pl.DataFrame({
        "d_allele": np.repeat(d_alleles, n5 * n3),
        "ndel5": np.tile(np.repeat(ndel5, n3), len(d_alleles)).astype(np.int16),
        "ndel3": np.tile(ndel3, n5 * len(d_alleles)).astype(np.int16),
        "p": np.tile(p, len(d_alleles)).astype(np.float64),
    })


def _placeholder_ins(ins_max: int) -> pl.DataFrame:
    length = np.arange(0, ins_max + 1)
    return pl.DataFrame({"length": length.astype(np.int16), "p": _decay_p(length.astype(float), 5.0)})


def _uniform_dinucl() -> pl.DataFrame:
    frm, to = np.repeat(np.arange(4), 4), np.tile(np.arange(4), 4)
    return pl.DataFrame({"from_nt": frm.astype(np.uint8), "to_nt": to.astype(np.uint8),
                         "p": np.full(16, 0.25)})


def from_arda(locus: str, organism: str = "human", *,
              palindrome_max: dict[str, int] | None = None, ins_max: int = 40) -> Model:
    """Build a recombination :class:`Model` whose gene set + germline come from **arda**.

    Genomic frames (names, CDR3-region germline, palindrome-extended ``cut_segment``,
    full germline + anchor for stitching, functionality) are sourced from arda via
    :mod:`vdjtools.model.reference` — so generated sequences carry arda's IMGT allele names
    and germline. The marginal ``tables`` are **placeholders** (uniform gene usage; small-trim /
    small-insert biased deletions/insertions with wide support) meant to be refit by
    :func:`vdjtools.model.infer.infer_native`; their ``ndel`` / ``length`` support ranges bound
    what EM can learn, so they are sized to each segment's full cut-segment length and ``ins_max``.

    Args:
        locus: e.g. ``"TRB"``, ``"IGH"``.
        organism: e.g. ``"human"``, ``"mouse"``.
        palindrome_max: Max palindromic nt per trimmable end (default matches OLGA human).
        ins_max: Maximum N-region insertion length in the placeholder tables.

    Returns:
        A validated :class:`Model` (VDJ if arda has D germline for the locus, else VJ).

    Raises:
        ImportError: If arda (the ``[model]`` extra) is not installed.
        ValueError: If arda has no germline for ``locus`` / ``organism``.
    """
    from . import reference as ref

    pal = palindrome_max or _DEFAULT_PALINDROME_MAX
    gl = ref.load_germline(locus, organism)
    full = ref.arda_full_germline(locus, organism)
    v_gl = gl.filter(pl.col("segment") == "V")
    j_gl = gl.filter(pl.col("segment") == "J")
    d_gl = gl.filter(pl.col("segment") == "D")
    chain_type = "VDJ" if d_gl.height else "VJ"

    genes_v = _arda_vj_genomic(v_gl, full, "v", pal["v_3"])
    genes_j = _arda_vj_genomic(j_gl, full, "j", pal["j_5"])
    v_alleles = genes_v["v_allele"].to_list()
    j_alleles = genes_j["j_allele"].to_list()
    v_cut = max((len(x) for x in genes_v["cut_segment"]), default=1)
    j_cut = max((len(x) for x in genes_j["cut_segment"]), default=1)

    tables: dict[str, pl.DataFrame] = {
        "v_choice": _uniform_choice(v_alleles, genes_v["functional"].to_list(), "v_allele"),
        "v_3_del": _placeholder_deletion(v_alleles, -pal["v_3"], v_cut, "v_allele"),
        "j_5_del": _placeholder_deletion(j_alleles, -pal["j_5"], j_cut, "j_allele"),
    }

    if chain_type == "VDJ":
        genes_d = _arda_d_genomic(d_gl, pal["d_5"], pal["d_3"])
        d_alleles = genes_d["d_allele"].to_list()
        d_cut = max((len(x) for x in genes_d["cut_segment"]), default=1)
        tables["j_choice"] = _uniform_choice(j_alleles, genes_j["functional"].to_list(), "j_allele")
        tables["d_gene"] = _uniform_cond_choice(j_alleles, d_alleles, genes_d["functional"].to_list(),
                                                "j_allele", "d_allele")
        tables["n_d"] = pl.DataFrame({"n_d": np.array([1], np.uint8), "p": np.array([1.0])})
        tables["d_del"] = _placeholder_d_del(d_alleles, -pal["d_5"], d_cut, -pal["d_3"], d_cut)
        tables["vd_ins"] = _placeholder_ins(ins_max)
        tables["dj_ins"] = _placeholder_ins(ins_max)
        tables["vd_dinucl"] = _uniform_dinucl()
        tables["dj_dinucl"] = _uniform_dinucl()
        events = {
            "v_choice": Event("v_choice", EventKind.GENE_CHOICE),
            "j_choice": Event("j_choice", EventKind.GENE_CHOICE),
            "d_gene": Event("d_gene", EventKind.GENE_CHOICE, ("j_choice",)),
            "n_d": Event("n_d", EventKind.N_D),
            "v_3_del": Event("v_3_del", EventKind.DELETION, ("v_choice",)),
            "j_5_del": Event("j_5_del", EventKind.DELETION, ("j_choice",)),
            "d_del": Event("d_del", EventKind.DELETION_2D, ("d_gene",)),
            "vd_ins": Event("vd_ins", EventKind.INS_LENGTH),
            "dj_ins": Event("dj_ins", EventKind.INS_LENGTH),
            "vd_dinucl": Event("vd_dinucl", EventKind.DINUCLEOTIDE),
            "dj_dinucl": Event("dj_dinucl", EventKind.DINUCLEOTIDE),
        }
        genomic = {"genes_v": genes_v, "genes_d": genes_d, "genes_j": genes_j}
    else:  # VJ
        tables["j_choice"] = _uniform_cond_choice(v_alleles, j_alleles, genes_j["functional"].to_list(),
                                                  "v_allele", "j_allele")
        tables["vj_ins"] = _placeholder_ins(ins_max)
        tables["vj_dinucl"] = _uniform_dinucl()
        events = {
            "v_choice": Event("v_choice", EventKind.GENE_CHOICE),
            "j_choice": Event("j_choice", EventKind.GENE_CHOICE, ("v_choice",)),
            "v_3_del": Event("v_3_del", EventKind.DELETION, ("v_choice",)),
            "j_5_del": Event("j_5_del", EventKind.DELETION, ("j_choice",)),
            "vj_ins": Event("vj_ins", EventKind.INS_LENGTH),
            "vj_dinucl": Event("vj_dinucl", EventKind.DINUCLEOTIDE),
        }
        genomic = {"genes_v": genes_v, "genes_j": genes_j}

    manifest = Manifest(locus=locus, organism=organism, chain_type=chain_type, events=events,
                        palindrome_max={k: int(v) for k, v in pal.items()}, source=f"arda:{locus}")
    return Model(manifest=manifest, tables=tables, genomic=genomic).validate()


def save_model(model: Model, path: str | Path) -> None:
    """Write a model to ``path/`` as ``manifest.json`` + one parquet per event/germline table."""
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    (out / "manifest.json").write_text(model.manifest.to_json())
    for name, df in model.tables.items():
        df.write_parquet(out / f"{name}.parquet")
    for name, df in model.genomic.items():
        df.write_parquet(out / f"{name}.parquet")


def load_model(path: str | Path) -> Model:
    """Load a native model directory written by :func:`save_model`."""
    src = Path(path)
    manifest = Manifest.from_json((src / "manifest.json").read_text())
    tables = {name: pl.read_parquet(src / f"{name}.parquet") for name in manifest.events}
    genomic_names = ["genes_v", "genes_j"] + (["genes_d"] if manifest.chain_type == "VDJ" else [])
    genomic = {name: pl.read_parquet(src / f"{name}.parquet") for name in genomic_names}
    return Model(manifest=manifest, tables=tables, genomic=genomic)
