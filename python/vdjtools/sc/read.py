"""10x / AIRR-Cell ingestion → a long Rearrangement frame keyed on ``cell_id``.

The load-bearing single-cell representation here is deliberately **flat**: one row
per productive contig, canonical AIRR Rearrangement columns, plus a ``cell_id`` (the
cell barcode) tying a cell's chains together and ``umi_count`` / ``clone_id`` carried
alongside. Everything downstream (pairing, QC, clustering) consumes that frame; the
AIRR *Data File* export (:func:`write_airr_cell`) is a secondary interchange layer.

10x CellRanger ``all_contig_annotations.csv`` names the junction column ``cdr3`` (with
the conserved Cys/Phe-Trp anchors included) — byte-identical to our ``cdr3_aa`` /
AIRR ``junction_aa`` convention (see ``io.schema``) — and ``cdr3_nt`` for the
nucleotide junction, so they map straight across.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from ..io.schema import (
    C_CALL,
    CDR3_AA,
    CDR3_NT,
    COUNT,
    D_CALL,
    J_CALL,
    LOCUS,
    V_CALL,
)

#: Receptor loci 10x can call; ``Multi`` / ``None`` / anything else is dropped.
VALID_LOCI: tuple[str, ...] = ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL")

CELL_ID = "cell_id"
SEQUENCE_ID = "sequence_id"
UMI_COUNT = "umi_count"
CLONE_ID = "clone_id"

#: Canonical single-cell long-frame columns, in order.
SC_COLUMNS: list[str] = [
    CELL_ID, SEQUENCE_ID, LOCUS,
    V_CALL, D_CALL, J_CALL, C_CALL,
    CDR3_AA, CDR3_NT, COUNT, UMI_COUNT, CLONE_ID,
]

_TRUTHY = ("1", "true", "t", "yes", "y")


def _truthy_expr(col: str) -> pl.Expr:
    """Boolean expression: is ``col`` one of the truthy tokens (case-insensitive)?"""
    return (pl.col(col).cast(pl.Utf8).str.strip_chars().str.to_lowercase()
            .is_in(list(_TRUTHY)))


def _pick(columns, *candidates: str) -> str | None:
    """Return the first of ``candidates`` present in ``columns`` (else ``None``)."""
    have = set(columns)
    for cand in candidates:
        if cand in have:
            return cand
    return None


def _read_csv_str(path: str | Path) -> pl.DataFrame:
    """Read a (possibly gzipped) 10x CSV with every column as a string."""
    return pl.read_csv(
        Path(path), separator=",", infer_schema_length=0,
        null_values=["", "NA", "None"], truncate_ragged_lines=True,
    )


def read_10x(
    all_contig: str | Path,
    consensus: str | Path | None = None,
    *,
    require_cell: bool = True,
    require_high_conf: bool = True,
) -> pl.DataFrame:
    """Read 10x ``all_contig_annotations.csv[.gz]`` into the canonical sc long frame.

    Keeps only productive, cell-associated, high-confidence contigs on a real receptor
    locus (:data:`VALID_LOCI`) with a resolved consensus id, one row per contig. When a
    ``consensus_annotations`` file is supplied, the per-cell contig's V/D/J calls are
    replaced by the matched consensus calls (joined on
    ``(raw_clonotype_id, raw_consensus_id) == (clonotype_id, consensus_id)``).

    Args:
        all_contig: Path to ``all_contig_annotations.csv`` (``.gz`` accepted). Columns
            follow CellRanger VDJ (``barcode, is_cell, contig_id, high_confidence,
            chain, v_gene, d_gene, j_gene, c_gene, productive, cdr3, cdr3_nt, reads,
            umis, raw_clonotype_id, raw_consensus_id``). ``*_call`` spellings are also
            accepted in place of ``*_gene``.
        consensus: Optional ``consensus_annotations.csv[.gz]`` to source consensus
            V/D/J calls from; if ``None`` the contig's own calls are used.
        require_cell: Drop contigs whose ``is_cell`` is not truthy (default ``True``).
        require_high_conf: Drop contigs whose ``high_confidence`` is not truthy
            (default ``True``).

    Returns:
        A ``pl.DataFrame`` with columns ``cell_id, sequence_id, locus, v_call, d_call,
        j_call, c_call, cdr3_aa, cdr3_nt, duplicate_count, umi_count, clone_id`` — one
        row per surviving productive contig.

    Raises:
        ValueError: If a required column (``barcode``, ``contig_id``, ``chain``,
            ``cdr3``, ``raw_consensus_id``) is missing from ``all_contig``.
    """
    df = _read_csv_str(all_contig)
    cols = df.columns
    for required in ("barcode", "contig_id", "chain", "cdr3", "raw_consensus_id"):
        if required not in cols:
            raise ValueError(f"read_10x: {all_contig!r} missing required column {required!r}")

    # Row filters: cell / high-confidence / productive / valid locus / has consensus id.
    if require_cell and "is_cell" in cols:
        df = df.filter(_truthy_expr("is_cell"))
    if require_high_conf and "high_confidence" in cols:
        df = df.filter(_truthy_expr("high_confidence"))
    if "productive" in cols:
        df = df.filter(_truthy_expr("productive"))
    df = df.filter(
        pl.col("chain").is_in(list(VALID_LOCI))
        & pl.col("raw_consensus_id").is_not_null()
    )

    v_col = _pick(cols, "v_gene", "v_call")
    d_col = _pick(cols, "d_gene", "d_call")
    j_col = _pick(cols, "j_gene", "j_call")
    c_col = _pick(cols, "c_gene", "c_call")
    reads_col = _pick(cols, "reads", "duplicate_count")
    umis_col = _pick(cols, "umis", "umi_count")
    clone_col = _pick(cols, "raw_clonotype_id", "clonotype_id")

    def _seg(col: str | None) -> pl.Expr:
        return pl.col(col) if col else pl.lit(None, dtype=pl.Utf8)

    def _count(col: str | None) -> pl.Expr:
        if col:
            return pl.col(col).cast(pl.Int64, strict=False)
        return pl.lit(None, dtype=pl.Int64)

    out = df.select(
        pl.col("barcode").alias(CELL_ID),
        pl.col("contig_id").alias(SEQUENCE_ID),
        pl.col("chain").alias(LOCUS),
        _seg(v_col).alias(V_CALL),
        _seg(d_col).alias(D_CALL),
        _seg(j_col).alias(J_CALL),
        _seg(c_col).alias(C_CALL),
        pl.col("cdr3").alias(CDR3_AA),
        (pl.col("cdr3_nt") if "cdr3_nt" in cols else pl.lit(None, dtype=pl.Utf8)).alias(CDR3_NT),
        _count(reads_col).alias(COUNT),
        _count(umis_col).alias(UMI_COUNT),
        (pl.col(clone_col) if clone_col else pl.lit(None, dtype=pl.Utf8)).alias(CLONE_ID),
        (pl.col("raw_consensus_id")).alias("_consensus_id"),
    )

    if consensus is not None:
        out = _join_consensus(out, consensus)

    return out.select(SC_COLUMNS)


def _join_consensus(contigs: pl.DataFrame, consensus: str | Path) -> pl.DataFrame:
    """Overwrite V/D/J calls from the matched consensus record (kept null-safe)."""
    cons = _read_csv_str(consensus)
    ccols = cons.columns
    cid = _pick(ccols, "clonotype_id")
    sid = _pick(ccols, "consensus_id")
    if cid is None or sid is None:
        # Nothing to join on — leave contig calls untouched.
        return contigs

    v_col = _pick(ccols, "v_gene", "v_call")
    d_col = _pick(ccols, "d_gene", "d_call")
    j_col = _pick(ccols, "j_gene", "j_call")
    keep = cons.select(
        pl.col(cid).alias(CLONE_ID),
        pl.col(sid).alias("_consensus_id"),
        (pl.col(v_col) if v_col else pl.lit(None, dtype=pl.Utf8)).alias("_cons_v"),
        (pl.col(d_col) if d_col else pl.lit(None, dtype=pl.Utf8)).alias("_cons_d"),
        (pl.col(j_col) if j_col else pl.lit(None, dtype=pl.Utf8)).alias("_cons_j"),
    ).unique(subset=[CLONE_ID, "_consensus_id"])

    joined = contigs.join(keep, on=[CLONE_ID, "_consensus_id"], how="left")
    return joined.with_columns(
        pl.coalesce("_cons_v", V_CALL).alias(V_CALL),
        pl.coalesce("_cons_d", D_CALL).alias(D_CALL),
        pl.coalesce("_cons_j", J_CALL).alias(J_CALL),
    ).drop("_cons_v", "_cons_d", "_cons_j")


def read_airr_cell(path: str | Path) -> pl.DataFrame:
    """Read an AIRR Rearrangement TSV carrying a ``cell_id`` column into the sc frame.

    Args:
        path: Path to an AIRR Rearrangement TSV (``.gz`` accepted) with at least a
            ``cell_id`` column plus the usual AIRR fields.

    Returns:
        A ``pl.DataFrame`` in the canonical sc long-frame layout (:data:`SC_COLUMNS`);
        columns absent from the file are filled with nulls. ``junction_aa`` /
        ``junction`` are accepted as sources for ``cdr3_aa`` / ``cdr3_nt``.

    Raises:
        ValueError: If the file has no ``cell_id`` column.
    """
    df = pl.read_csv(Path(path), separator="\t", infer_schema_length=0,
                     null_values=["", "NA", "None"], truncate_ragged_lines=True)
    cols = df.columns
    if CELL_ID not in cols:
        raise ValueError(f"read_airr_cell: {path!r} has no 'cell_id' column")

    def _str(name, *alts) -> pl.Expr:
        src = _pick(cols, name, *alts)
        return (pl.col(src) if src else pl.lit(None, dtype=pl.Utf8)).alias(name)

    def _int(name, *alts) -> pl.Expr:
        src = _pick(cols, name, *alts)
        return (pl.col(src).cast(pl.Int64, strict=False) if src
                else pl.lit(None, dtype=pl.Int64)).alias(name)

    if LOCUS in cols:
        locus_expr = pl.col(LOCUS)
    elif V_CALL in cols:
        locus_expr = pl.col(V_CALL).str.slice(0, 3)  # derive locus from v_call prefix
    else:
        locus_expr = pl.lit(None, dtype=pl.Utf8)
    return df.select(
        _str(CELL_ID),
        _str(SEQUENCE_ID),
        locus_expr.alias(LOCUS),
        _str(V_CALL), _str(D_CALL), _str(J_CALL), _str(C_CALL),
        _str(CDR3_AA, "junction_aa"), _str(CDR3_NT, "junction"),
        _int(COUNT, "reads"), _int(UMI_COUNT, "umis"),
        _str(CLONE_ID, "raw_clonotype_id", "clonotype_id"),
    )


def _receptor_hash(dom1_aa: str, dom2_aa: str) -> str:
    """AIRR ``receptor_hash`` = ``sha256(upper(dom1) + upper(dom2))`` hex, lowercased."""
    import hashlib

    payload = (str(dom1_aa or "").upper() + str(dom2_aa or "").upper()).encode()
    return hashlib.sha256(payload).hexdigest()


def write_airr_cell(
    rearr: pl.DataFrame,
    cells_out: str | Path,
    *,
    receptors: bool = True,
    repertoire_id: str = "",
) -> Path:
    """Emit an AIRR Data File (YAML) with a ``Cell`` array (and optional ``Receptor``).

    Builds one ``Cell`` per ``cell_id`` (linking its ``sequence_id`` receptors) and,
    when ``receptors`` is set, one ``Receptor`` per paired heavy/light chain within a
    cell. The ``receptor_hash`` is ``sha256`` of the two upper-cased domain sequences.

    .. note::
        The AIRR spec's ``receptor_variable_domain_{1,2}_aa`` is the **full mature
        V-domain** amino-acid sequence. 10x contigs only expose the junction, so this
        field is populated with the **junction** (``cdr3_aa``) and that limitation is
        recorded in the file's ``Info`` block. Downstream code should treat these as
        junction-level, not full-domain, sequences.

    Args:
        rearr: A single-cell long frame (:data:`SC_COLUMNS`), typically the paired
            output — one ``cell_id`` may carry several chains.
        cells_out: Destination path for the AIRR Data File (``.yaml`` / ``.json``).
        receptors: Emit the ``Receptor`` list pairing heavy (β/heavy) and light
            (α/light) chains per cell (default ``True``).
        repertoire_id: Value stamped into each ``Cell.repertoire_id`` (default empty).

    Returns:
        The ``Path`` written.

    Raises:
        ImportError: If PyYAML is not installed (see the ``sc`` extra).
    """
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover - exercised only without pyyaml
        raise ImportError(
            "PyYAML is required for write_airr_cell; install the extra with "
            "`pip install 'vdjtools[sc]'` (or `pip install pyyaml`)."
        ) from exc

    heavy = {"TRB", "TRD", "IGH"}
    cell_rows: list[dict] = []
    receptor_rows: list[dict] = []

    for cell_id, grp in rearr.group_by(CELL_ID, maintain_order=True):
        cid = cell_id[0] if isinstance(cell_id, tuple) else cell_id
        contigs = grp.to_dicts()
        receptor_ids: list[str] = []

        if receptors:
            heavies = [c for c in contigs if (c.get(LOCUS) or "") in heavy]
            lights = [c for c in contigs if (c.get(LOCUS) or "") not in heavy]
            for h in heavies:
                for lt in lights:
                    dom1 = h.get(CDR3_AA) or ""
                    dom2 = lt.get(CDR3_AA) or ""
                    rid = f"{cid}:{h.get(SEQUENCE_ID)}:{lt.get(SEQUENCE_ID)}"
                    receptor_ids.append(rid)
                    rtype = "TCR" if (h.get(LOCUS) or "").startswith("TR") else "BCR"
                    receptor_rows.append({
                        "receptor_id": rid,
                        "receptor_hash": _receptor_hash(dom1, dom2),
                        "receptor_type": rtype,
                        "receptor_variable_domain_1_aa": dom1,
                        "receptor_variable_domain_1_locus": h.get(LOCUS),
                        "receptor_variable_domain_2_aa": dom2,
                        "receptor_variable_domain_2_locus": lt.get(LOCUS),
                    })
        if not receptor_ids:
            receptor_ids = [str(c.get(SEQUENCE_ID)) for c in contigs]

        cell_rows.append({
            "cell_id": str(cid),
            "repertoire_id": repertoire_id,
            "virtual_pairing": False,
            "receptors": receptor_ids,
        })

    doc = {
        "Info": {
            "title": "vdjtools single-cell export",
            "note": (
                "receptor_variable_domain_*_aa hold the CDR3/junction (anchors "
                "included), NOT the full mature V-domain: 10x contigs expose only the "
                "junction. Treat these as junction-level sequences."
            ),
        },
        "Cell": cell_rows,
    }
    if receptors:
        doc["Receptor"] = receptor_rows

    out = Path(cells_out)
    out.write_text(yaml.safe_dump(doc, sort_keys=False))
    return out
