"""10x / AIRR-Cell ingestion → a long Rearrangement frame keyed on ``cell_id``.

The load-bearing single-cell representation here is deliberately **flat**: one row
per productive contig, canonical AIRR Rearrangement columns, plus a ``cell_id`` (the
cell barcode) tying a cell's chains together and ``umi_count`` / ``clone_id`` carried
alongside. Everything downstream (pairing, QC, clustering) consumes that frame; the
AIRR *Data File* export (:func:`write_airr_cell`) is a secondary interchange layer.

10x CellRanger ``all_contig_annotations.csv`` names the junction column ``cdr3`` (with
the conserved Cys/Phe-Trp anchors included) — content-identical to our ``junction_aa`` /
AIRR ``junction_aa`` convention (see ``io.schema``) — and ``cdr3_nt`` for the
nucleotide junction; the reader maps those straight onto ``junction_aa`` / ``junction_nt``.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from ..io.schema import (
    C_CALL,
    JUNCTION_AA,
    JUNCTION_NT,
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
PRODUCTIVE = "productive"

#: Canonical single-cell long-frame columns, in order.
#:
#: ``productive`` is a mandatory AIRR Rearrangement field and is carried so
#: :func:`vdjtools.sc.to_airr` can emit a schema-valid table; the readers here keep only
#: productive contigs, so it is ``True`` wherever the source stated productivity at all.
SC_COLUMNS: list[str] = [
    CELL_ID, SEQUENCE_ID, LOCUS,
    V_CALL, D_CALL, J_CALL, C_CALL,
    JUNCTION_AA, JUNCTION_NT, COUNT, UMI_COUNT, CLONE_ID, PRODUCTIVE,
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
    """Read a 10x contig-annotation CSV into the canonical sc long frame.

    Accepts both ``all_contig_annotations.csv[.gz]`` and
    ``filtered_contig_annotations.csv[.gz]`` — the two share one writer in CellRanger and
    so one column layout; ``filtered_`` is simply pre-restricted to
    ``is_cell && high_confidence``, which this reader applies anyway.

    Keeps only productive, cell-associated, high-confidence contigs on a real receptor
    locus (:data:`VALID_LOCI`) with a resolved consensus id, one row per contig. When a
    ``consensus_annotations`` file is supplied, the per-cell contig's V/D/J calls are
    replaced by the matched consensus calls (joined on
    ``(raw_clonotype_id, raw_consensus_id) == (clonotype_id, consensus_id)``).

    Column drift across CellRanger versions is tolerated: the ``fwr*``/``cdr1``/``cdr2``
    region columns exist only from CR6, ``exact_subclonotype_id`` from CR4+, and ``sample``
    only in ``cellranger multi`` output — none is required here. ``raw_consensus_id`` is
    used when present and skipped when absent.

    NOTE: For CellRanger 4.0+, prefer :func:`read_airr_cell` on ``airr_rearrangement.tsv``
    — it is the same file the downstream tools read, so ingestion and interop agree.

    Args:
        all_contig: Path to ``all_contig_annotations.csv`` or
            ``filtered_contig_annotations.csv`` (``.gz`` accepted). Columns follow
            CellRanger VDJ (``barcode, is_cell, contig_id, high_confidence,
            chain, v_gene, d_gene, j_gene, c_gene, productive, cdr3, cdr3_nt, reads,
            umis, raw_clonotype_id, raw_consensus_id``). ``*_call`` spellings are also
            accepted in place of ``*_gene``.
        consensus: Optional ``consensus_annotations.csv[.gz]`` to source consensus
            V/D/J calls from; if ``None`` the contig's own calls are used.
        require_cell: Drop contigs whose ``is_cell`` is not truthy (default ``True``).
        require_high_conf: Drop contigs whose ``high_confidence`` is not truthy
            (default ``True``).

    Returns:
        A ``pl.DataFrame`` in the canonical sc long-frame layout (:data:`SC_COLUMNS`) — one
        row per surviving productive contig.

    Raises:
        ValueError: If a required column (``barcode``, ``contig_id``, ``chain``,
            ``cdr3``) is missing from ``all_contig``.
    """
    df = _read_csv_str(all_contig)
    cols = df.columns
    for required in ("barcode", "contig_id", "chain", "cdr3"):
        if required not in cols:
            raise ValueError(f"read_10x: {all_contig!r} missing required column {required!r}")

    # Row filters: cell / high-confidence / productive / valid locus / has consensus id.
    if require_cell and "is_cell" in cols:
        df = df.filter(_truthy_expr("is_cell"))
    if require_high_conf and "high_confidence" in cols:
        df = df.filter(_truthy_expr("high_confidence"))
    if "productive" in cols:
        df = df.filter(_truthy_expr("productive"))
    df = df.filter(pl.col("chain").is_in(list(VALID_LOCI)))
    # `raw_consensus_id` is absent from some CellRanger versions/outputs; only gate on it
    # when it is actually there (an unassigned contig has no consensus call to trust).
    if "raw_consensus_id" in cols:
        df = df.filter(pl.col("raw_consensus_id").is_not_null())

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
        pl.col("cdr3").alias(JUNCTION_AA),
        (pl.col("cdr3_nt") if "cdr3_nt" in cols else pl.lit(None, dtype=pl.Utf8)).alias(JUNCTION_NT),
        _count(reads_col).alias(COUNT),
        _count(umis_col).alias(UMI_COUNT),
        (pl.col(clone_col) if clone_col else pl.lit(None, dtype=pl.Utf8)).alias(CLONE_ID),
        # Non-productive contigs were filtered out above, so this is True by construction
        # (null only when the file never stated productivity).
        pl.lit(True if "productive" in cols else None, dtype=pl.Boolean).alias(PRODUCTIVE),
        (pl.col("raw_consensus_id") if "raw_consensus_id" in cols
         else pl.lit(None, dtype=pl.Utf8)).alias("_consensus_id"),
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

    This is also the reader for **CellRanger's ``airr_rearrangement.tsv``** (emitted since
    Cell Ranger 4.0, and under ``per_sample_outs/<sample>/vdj_t/`` for ``cellranger multi``)
    and for **arda's** barcoded output — both are plain AIRR Rearrangement tables with a
    ``cell_id`` column. Prefer it over :func:`read_10x`: it is the same file scirpy,
    dandelion and scRepertoire read, so the ingestion and interop paths agree by
    construction.

    Args:
        path: Path to an AIRR Rearrangement TSV (``.gz`` accepted) with at least a
            ``cell_id`` column plus the usual AIRR fields.

    Returns:
        A ``pl.DataFrame`` in the canonical sc long-frame layout (:data:`SC_COLUMNS`);
        columns absent from the file are filled with nulls. ``junction_aa`` /
        ``junction`` are accepted as sources for ``junction_aa`` / ``junction_nt``.

    Raises:
        ValueError: If the file has no ``cell_id`` column.
    """
    df = pl.read_csv(Path(path), separator="\t", infer_schema_length=0,
                     null_values=["", "NA", "None"], truncate_ragged_lines=True)
    if CELL_ID not in df.columns:
        raise ValueError(f"read_airr_cell: {path!r} has no 'cell_id' column")
    return _airr_frame_to_sc(df)


def _airr_frame_to_sc(df: pl.DataFrame) -> pl.DataFrame:
    """Map an already-loaded AIRR Rearrangement frame onto :data:`SC_COLUMNS`.

    Shared by :func:`read_airr_cell` and :func:`read_arda_cells`, which differ only in how
    they get the frame off disk.
    """
    cols = df.columns

    def _str(name, *alts, alias=None) -> pl.Expr:
        src = _pick(cols, name, *alts)
        return (pl.col(src) if src else pl.lit(None, dtype=pl.Utf8)).alias(alias or name)

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
        # Prefer the junction (anchors INCLUDED) over IMGT cdr3_aa/cdr3_nt (excluded),
        # matching io/read.py and the canonical junction_aa=junction convention.
        _str("junction_aa", "cdr3_aa", alias=JUNCTION_AA),
        _str("junction_nt", "junction", "cdr3_nt", "cdr3", alias=JUNCTION_NT),
        _int(COUNT, "reads", "consensus_count"), _int(UMI_COUNT, "umis"),
        _str(CLONE_ID, "raw_clonotype_id", "clonotype_id"),
        (_truthy_expr(PRODUCTIVE) if PRODUCTIVE in cols
         else pl.lit(None, dtype=pl.Boolean)).alias(PRODUCTIVE),
    )


#: Cell-level QC columns lifted from arda's ``.chains.tsv``, namespaced so they cannot be
#: confused with :func:`vdjtools.sc.resolve_chains`' own verdict.
ARDA_STATUS = "arda_status"
ARDA_MOLECULES = "arda_molecules"


def read_arda_cells(prefix: str | Path, *, chains: bool = True) -> pl.DataFrame:
    """Read the output of arda's single-cell pipeline (``arda cells``).

    ``arda cells`` writes ``<prefix>.contigs.airr.tsv`` — an AIRR Rearrangement table with
    ``cell_id``, ``molecules`` and ``reads`` — plus a cell-level ``<prefix>.chains.tsv``
    carrying arda's own per-chain verdict (``status`` is one of ``primary``, ``secondary``,
    ``doublet_candidate``, ``extra``).

    NOTE: arda's ``status`` and :func:`~vdjtools.sc.resolve_chains`' verdict are
    **independent** calls on the same question. This reader surfaces arda's as
    ``arda_status`` rather than acting on it, so running ``resolve_chains`` afterwards does
    not silently discard the upstream judgement -- compare them, don't assume they agree.

    Reading goes through arda's own ``read_airr``, not a plain CSV read: arda has written
    two TSV dialects and that reader is what normalises them. It matters because
    ``junction_quality`` is Phred+33 and character 34 is a double quote.

    Args:
        prefix: The ``arda cells`` output prefix (``<prefix>.contigs.airr.tsv`` is read), or
            a direct path to a ``*.airr.tsv`` file.
        chains: Join arda's per-chain ``status`` / ``molecules`` from ``<prefix>.chains.tsv``
            when that file exists (default ``True``).

    Returns:
        A ``pl.DataFrame`` in the canonical layout (:data:`SC_COLUMNS`), plus
        ``arda_status`` / ``arda_molecules`` when the chains table was joined.

    Raises:
        FileNotFoundError: If no contigs AIRR table is found for ``prefix``.
    """
    from arda.annotate.airr_out import read_airr as _arda_read_airr

    p = Path(prefix)
    contigs = p if p.name.endswith(".airr.tsv") else Path(f"{p}.contigs.airr.tsv")
    if not contigs.exists():
        raise FileNotFoundError(
            f"no arda contigs table at {contigs} -- expected `arda cells -p {prefix}` output"
        )

    df = _arda_read_airr(contigs)
    if CELL_ID not in df.columns:
        raise ValueError(
            f"{contigs} has no 'cell_id' column; for arda's BULK output use "
            "vdjtools.io.read_arda, or re-run with `arda map --cell-from/--cell-regex`"
        )
    # arda names the per-contig molecule count `molecules`; the canonical frame wants it as
    # umi_count (one molecule == one UMI-tagged consensus).
    if "molecules" in df.columns and UMI_COUNT not in df.columns:
        df = df.rename({"molecules": UMI_COUNT})
    out = _airr_frame_to_sc(df)

    chains_tsv = Path(f"{p}.chains.tsv")
    if chains and chains_tsv.exists():
        out = _join_arda_chains(out, chains_tsv)
    return out


def _join_arda_chains(cells: pl.DataFrame, chains_tsv: Path) -> pl.DataFrame:
    """Attach arda's per-chain ``status`` / ``molecules`` on (cell_id, locus, junction_aa)."""
    ch = pl.read_csv(chains_tsv, separator="\t", infer_schema_length=0,
                     null_values=["", "NA", "None"], truncate_ragged_lines=True)
    key = [CELL_ID, LOCUS, JUNCTION_AA]
    if not set(key) <= set(ch.columns) or "status" not in ch.columns:
        return cells
    keep = ch.select(
        *key,
        pl.col("status").alias(ARDA_STATUS),
        (pl.col("molecules").cast(pl.Int64, strict=False) if "molecules" in ch.columns
         else pl.lit(None, dtype=pl.Int64)).alias(ARDA_MOLECULES),
    ).unique(subset=key, keep="first")
    return cells.join(keep, on=key, how="left")


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
        field is populated with the **junction** (``junction_aa``) and that limitation is
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
                    dom1 = h.get(JUNCTION_AA) or ""
                    dom2 = lt.get(JUNCTION_AA) or ""
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
