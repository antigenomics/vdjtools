"""Converters for legacy repertoire input formats → the canonical clonotype frame.

Reimplements the battle-tested legacy vdjtools format parsers (their exact column mappings,
gene-name normalisation and CDR3 handling were lifted verbatim from the Groovy
``com.antigenomics.vdjtools.io.parser`` classes) so third-party tool output can be read
straight into the canonical AIRR-junction frame of :mod:`vdjtools.io.schema`:

* **MiXcr** (v1/2 *and* v3/4 header dialects, incl. the C-gene / isotype hit) — :func:`read_mixcr`
* **MiGec** — :func:`read_migec`
* **Adaptive immunoSEQ** v1 and v2 — :func:`read_immunoseq`
* **IMGT/HighV-QUEST** — :func:`read_imgt`
* **Vidjil** (``.vidjil`` JSON) — :func:`read_vidjil`
* **RTCR** — :func:`read_rtcr`
* **TRUST4** (``*_report.tsv``) — :func:`read_trust4`
* **arda** (AIRR annotation output) — :func:`read_arda`

Every reader returns the canonical frame: V/D/J IMGT calls, ``junction_nt`` / ``junction_aa``
(the AIRR junction — conserved Cys104 … Phe/Trp118 anchors **included**, matching the legacy
vdjtools ``cdr3nt``/``cdr3aa``), ``duplicate_count`` and recomputed ``frequency``. Per-read
formats are collapsed to unique clonotypes with summed counts.
"""
from __future__ import annotations

import csv
import functools
import gzip
import json
import os
import re
from importlib import resources
from pathlib import Path

import polars as pl

from . import schema
from .read import _read_tsv, read_airr
from .schema import (
    C_CALL,
    COUNT,
    D_CALL,
    J_CALL,
    JUNCTION_AA,
    JUNCTION_NT,
    V_CALL,
)


def _to_int(*cells) -> int:
    """First cell parsing to a **positive** integer count (double-then-truncate); none → 0.

    Skips non-numeric *and* non-positive cells so the legacy count fallback fires correctly —
    e.g. Adaptive immunoSEQ v1 ``templates`` is often ``"null"`` **or** ``"0"`` and the real
    count lives in ``reads``. A genuinely zero count is dropped downstream by :func:`_finalize`.
    """
    for x in cells:
        try:
            v = int(float(x))
        except (TypeError, ValueError):
            continue
        if v > 0:
            return v
    return 0


def _read_text(path: str | os.PathLike) -> str:
    """Read a text file, transparently decompressing gzip (``.gz`` or magic bytes)."""
    data = Path(path).read_bytes()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode()

# --- codon table + legacy-faithful translation ------------------------------------------

_CODON2AA = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L",
    "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*",
    "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C",
    "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G",
}


def _codon2aa(codon: str) -> str:
    return _CODON2AA.get(codon.upper(), "X")


def translate(seq: str) -> str:
    """Translate a (possibly out-of-frame) CDR3 nt sequence, V→J bidirectionally.

    Port of the legacy ``CommonUtil.translate``: in-frame sequences are a plain codon walk
    (stop → ``*``); an out-of-frame sequence is padded in the middle with ``?`` and translated
    inward from both ends, leaving the untranslatable middle codon(s) lower-cased (later
    collapsed to ``_`` by :func:`to_unified_cdr3aa`).
    """
    if not seq:
        return ""
    oof = len(seq) % 3
    if oof:
        mid = len(seq) // 2
        seq = seq[:mid] + "?" * (3 - oof) + seq[mid:]

    left_end = right_end = -1
    aa = ""
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if "?" in codon:
            left_end = i
            break
        aa += _codon2aa(codon)
    if oof == 0:
        return aa

    aa_right = ""
    for i in range(len(seq), 2, -3):
        codon = seq[i - 3:i]
        if "?" in codon:
            right_end = i
            break
        aa_right += _codon2aa(codon)
    return aa + seq[left_end:right_end].lower() + aa_right[::-1]


_OOF_RUN = re.compile(r"[atgc#~_?]+")


def to_unified_cdr3aa(seq: str | None) -> str | None:
    """Collapse each run of non-coding markers (lower-case nt / ``# ~ _ ?``) to a single ``_``."""
    if seq is None:
        return None
    return _OOF_RUN.sub("_", seq)


# --- gene-call normalisation ------------------------------------------------------------

def extract_vdj(field: str | None) -> str | None:
    """First tie, allele stripped, quotes/space trimmed → an IMGT gene call (or ``None``).

    Port of ``CommonUtil.extractVDJ``: ``"TRBV12-4,TRBV12-3" → "TRBV12-4"``,
    ``"TRBV13*00(401.5)" → "TRBV13"``. Empty / placeholder → ``None``.
    """
    if field is None:
        return None
    gene = field.split(",")[0].split("*")[0].replace('"', "").strip()
    return gene or None


_ZERO_PAD = re.compile(r"0([1-9])")


@functools.lru_cache(maxsize=1)
def _adaptive_map() -> dict[str, str]:
    """Adaptive token → IMGT gene, from the shipped CDR-validated table (cached).

    Built by ``appendix/build_adaptive_imgt_map.py``. The table's own ``evidence`` column
    carries, per row, how that token was resolved -- CDR-sequence agreement, a family
    fallback, or a validated manual call -- so the rationale travels with the data.
    """
    txt = resources.files("vdjtools.resources").joinpath("adaptive_imgt_map.tsv").read_text()
    return {r["adaptive_token"]: r["imgt_gene"]
            for r in csv.DictReader(txt.splitlines(), delimiter="\t")
            if r["imgt_gene"] and r["status"] != "unresolved_dropped"}


def _adaptive_to_imgt(field: str | None) -> str | None:
    """Adaptive immunoSEQ → IMGT gene name, via the shipped lookup table.

    ``"TCRBV29-01" → "TRBV29-1"``, ``"TCRAJ39-01" → "TRAJ39"``,
    ``"TCRAV38-02" → "TRAV38-2/DV8"``, ``"TCRBV12-X" → "TRBV12-3"``; ``"unresolved" → None``.

    Whether an Adaptive token's trailing group is an IMGT subgroup or an allele is a per-family
    fact (``TCRAV01-01`` = ``TRAV1-1`` but ``TCRAV22-01`` = ``TRAV22``), so no regex can decide it:
    the legacy zero-strip of ``CommonUtil.extractVDJImmunoSeq`` emits a *non-existent* gene name
    for 100 of the 161 tokens seen in IMMREP25 + pairSEQ. Off-table tokens fall back to that
    legacy rewrite, so unknown input behaves exactly as before.
    """
    if field is None:
        return None
    # look the raw token up before any splitting — slash ties (``TCRBV03-01/03-02``) are keys.
    tok = field.split(",")[0].replace('"', "").strip().split("*")[0]
    hit = _adaptive_map().get(tok)
    if hit:
        return hit
    gene = extract_vdj(field)
    if gene is None or gene.lower() == "unresolved":
        return None
    return _ZERO_PAD.sub(r"\1", gene.replace("TCR", "TR"))


def _adaptive_call(gene: str | None, family: str | None, ties: str | None) -> str | None:
    """Resolve an Adaptive call gene→family→familyTies (each ``_adaptive_to_imgt``)."""
    return _adaptive_to_imgt(gene) or _adaptive_to_imgt(family) or _adaptive_to_imgt(ties)


# --- shared finalisation ----------------------------------------------------------------

def _lower_map(cols: list[str]) -> dict[str, str]:
    return {c.lower(): c for c in cols}


def _pick(lower: dict[str, str], *names: str) -> str | None:
    """Return the actual column name for the first case-insensitive match, else ``None``."""
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    return None


#: Stands in for a null inside the Adaptive dedup join. Joining on null keys does not match in
#: polars, and the triples genuinely contain nulls, so they are carried as a byte no gene name has.
_NULL_KEY = "\x00"


def _adaptive_col(raw: pl.DataFrame, gene, family, ties) -> pl.Series:
    """:func:`_adaptive_call` over a column, once per DISTINCT ``(gene, family, ties)`` triple.

    The resolution is a shipped-table lookup with a legacy-rewrite fallback, so it is not one
    expression -- but it is a pure function of the triple, and an export of tens of thousands of
    rearrangements carries a few dozen distinct triples. Deduplication, not a cache: the table is
    built and dropped inside this call, and the mapping back onto the column is a join.
    """
    def key(col):
        e = pl.col(col).cast(pl.Utf8) if col else pl.lit(None, pl.Utf8)
        return e.fill_null(_NULL_KEY)

    names = ["_g", "_f", "_t"]
    keys = raw.select([key(c).alias(n) for c, n in zip((gene, family, ties), names)])
    uniq = keys.unique(maintain_order=True)
    resolved = [_adaptive_call(*(None if x == _NULL_KEY else x for x in row))
                for row in uniq.iter_rows()]
    table = uniq.with_columns(pl.Series("_call", resolved, dtype=pl.Utf8))
    return keys.join(table, on=names, how="left")["_call"]


def _vdj_expr(col: str | None) -> pl.Expr:
    """:func:`extract_vdj` as one polars expression -- a whole column in one pass, not per row.

    The same four steps in the same order (first comma-separated tie, allele stripped, quotes and
    surrounding space removed), and an empty result is null, matching the scalar form's ``None``.
    A ``None`` column name yields an all-null column, which is how an absent optional D/C is
    expressed without a branch at every call site.
    """
    if col is None:
        return pl.lit(None, pl.Utf8)
    gene = (pl.col(col).cast(pl.Utf8)
            .str.split(",").list.first()
            .str.split("*").list.first()
            .str.replace_all('"', "", literal=True)
            .str.strip_chars())
    return pl.when(gene.str.len_bytes() > 0).then(gene).otherwise(None)


def _upper_expr(col: str) -> pl.Expr:
    """A nucleotide column, upper-cased, with empty as null -- ``(x or "").upper() or None``."""
    up = pl.col(col).cast(pl.Utf8).fill_null("").str.to_uppercase()
    return pl.when(up.str.len_bytes() > 0).then(up).otherwise(None)


def _unified_expr(col: str) -> pl.Expr:
    """:func:`to_unified_cdr3aa` as an expression: each run of non-coding markers to one ``_``."""
    return pl.col(col).cast(pl.Utf8).str.replace_all(r"[atgc#~_?]+", "_")


def _translated(nt: pl.Series) -> pl.Series:
    """:func:`translate` then :func:`to_unified_cdr3aa` over a whole column, natively.

    The translation is a codon walk -- bidirectional on an out-of-frame junction, which polars
    cannot express -- and it was the last per-row Python in the readers: a 42,877-row immunoSEQ
    export walks ~40,000 of them, so deduplicating over distinct sequences saves almost nothing.
    ``_core.translate_junctions`` is the same function in C++, checked against this module's
    scalar pair on 4,009 sequences (every frame offset, with and without non-ACGT bases).

    A null or empty junction stays null, as the scalar form's ``if nt else None`` did.
    """
    from .._core import translate_junctions

    aa = pl.Series("aa", translate_junctions(nt.fill_null("").to_list()), dtype=pl.Utf8)
    blank = pl.Series("aa", [None] * nt.len(), dtype=pl.Utf8)
    return aa.zip_with(nt.is_not_null() & (nt.str.len_bytes() > 0), blank)


def _finalize(df: pl.DataFrame) -> pl.DataFrame:
    """Filter bad rows, collapse to unique clonotypes (summed counts), coerce to canonical.

    Takes the frame the reader built with expressions. It used to take a list of one dict per
    input row: on a 42,877-row immunoSEQ export that was 42,877 dict allocations plus a schema
    inference pass over them, before any of the work polars is for.
    """
    if df.height == 0:
        return schema.add_locus(schema.normalize(pl.DataFrame(schema={c: pl.Utf8 for c in
                                (V_CALL, D_CALL, J_CALL, JUNCTION_AA, JUNCTION_NT)})))
    df = df.with_columns(pl.col(COUNT).cast(pl.Int64, strict=False))
    keep = (
        pl.col(JUNCTION_NT).is_not_null() & (pl.col(JUNCTION_NT).str.len_bytes() > 0)
        & pl.col(JUNCTION_AA).is_not_null() & (pl.col(JUNCTION_AA).str.len_bytes() > 0)
        & pl.col(V_CALL).is_not_null() & pl.col(J_CALL).is_not_null()
        & pl.col(COUNT).is_not_null() & (pl.col(COUNT) > 0)
    )
    df = df.filter(keep)
    key = [V_CALL, J_CALL, JUNCTION_NT, JUNCTION_AA]
    reps = [pl.col(c).drop_nulls().first().alias(c) for c in (D_CALL, C_CALL) if c in df.columns]
    df = df.group_by(key, maintain_order=True).agg(pl.col(COUNT).sum(), *reps)
    return schema.add_locus(schema.normalize(df, recompute_freq=True))


# --- format readers ---------------------------------------------------------------------

def read_mixcr(path: str | os.PathLike, n_rows: int | None = None) -> pl.DataFrame:
    """Read a MiXcr ``exportClones`` table — legacy (v1/2) or current (v3/4) dialect.

    The two header dialects name every field differently, so both spellings are accepted
    for each column:

    ==============  ==========================  ================================================
    field           legacy (v1/2)               current (v3/4)
    ==============  ==========================  ================================================
    count           ``Clone count``             ``cloneCount`` / ``readCount`` / ``uniqueTagCountMolecule``
    V/D/J/C hits    ``All V hits`` …            ``allVHitsWithScore`` …
    CDR3 nt / aa    ``N. Seq. CDR3`` …          ``nSeqCDR3`` / ``aaSeqCDR3``
    ==============  ==========================  ================================================

    MiXcr's current field API exposes ``-readCount`` / ``-readFraction`` (there is no
    ``-cloneCount`` field any more), while the default preset still *labels* the column
    ``cloneCount`` — so a v4 export carries one spelling or the other depending on how it
    was produced, and both must parse.

    Count precedence is read-based first (``cloneCount`` → ``readCount``), with the UMI
    molecule count (``uniqueTagCountMolecule``) used only when no read count column is
    present. On a UMI library the molecule count is the less PCR-biased abundance, so
    prefer exporting it alone if that is the quantity you want.
    """
    raw = _read_tsv(path, n_rows=n_rows)
    lo = _lower_map(raw.columns)
    count_c = _pick(lo, "clone count", "clonecount", "readcount", "uniquetagcountmolecule")
    v_c = _pick(lo, "all v hits", "allvhitswithscore")
    d_c = _pick(lo, "all d hits", "alldhitswithscore")
    j_c = _pick(lo, "all j hits", "alljhitswithscore")
    c_c = _pick(lo, "all c hits", "allchitswithscore")  # C gene / BCR isotype (v1/2 & v3/4)
    nt_c = _pick(lo, "n. seq. cdr3", "nseqcdr3", "nseqimputedcdr3")
    aa_c = _pick(lo, "aa. seq. cdr3", "aaseqcdr3", "aaseqimputedcdr3")
    if not (count_c and v_c and j_c and nt_c and aa_c):
        raise ValueError(f"not a MiXcr table (need count / V,J hits / CDR3 nt+aa); have {raw.columns}")
    return _finalize(raw.select(
        _vdj_expr(v_c).alias(V_CALL), _vdj_expr(d_c).alias(D_CALL),
        _vdj_expr(j_c).alias(J_CALL), _vdj_expr(c_c).alias(C_CALL),
        _upper_expr(nt_c).alias(JUNCTION_NT),
        # MiXcr aa is milib-based -- kept verbatim (no unify)
        pl.col(aa_c).cast(pl.Utf8).alias(JUNCTION_AA),
        pl.col(count_c).cast(pl.Float64, strict=False).cast(pl.Int64).alias(COUNT),
    ))


def read_migec(path: str | os.PathLike, n_rows: int | None = None) -> pl.DataFrame:
    """Read a MiGEC ``CdrBlast`` clonotype table."""
    raw = _read_tsv(path, n_rows=n_rows)
    lo = _lower_map(raw.columns)
    count_c = _pick(lo, "count")
    nt_c = _pick(lo, "cdr3 nucleotide sequence")
    aa_c = _pick(lo, "cdr3 amino acid sequence")
    v_c = _pick(lo, "v segments")
    j_c = _pick(lo, "j segments")
    d_c = _pick(lo, "d segments")
    if not (count_c and nt_c and aa_c and v_c and j_c):
        raise ValueError(f"not a MiGEC table; have {raw.columns}")
    return _finalize(raw.select(
        _vdj_expr(v_c).alias(V_CALL), _vdj_expr(d_c).alias(D_CALL),
        _vdj_expr(j_c).alias(J_CALL),
        _upper_expr(nt_c).alias(JUNCTION_NT),
        _unified_expr(aa_c).alias(JUNCTION_AA),
        pl.col(count_c).cast(pl.Float64, strict=False).cast(pl.Int64).alias(COUNT),
    ))


def read_mitcr(path: str | os.PathLike, n_rows: int | None = None) -> pl.DataFrame:
    """Read a MiTCR / tcR (R package) clonotype table — the dot-separated dialect.

    Header is ``Read.count Read.proportion CDR3.nucleotide.sequence CDR3.amino.acid.sequence
    V.gene J.gene D.gene V.end J.start D5.end D3.end VD.insertions DJ.insertions
    Total.insertions``. Distinct from MiGEC's ``CDR3 nucleotide sequence`` / ``V segments``
    (spaces, not dots), so it needs its own picks. ``D.gene`` may carry an ambiguous call
    (``"TRBD1, TRBD2"``); :func:`extract_vdj` keeps the first.

    Args:
        path: Path to a MiTCR/tcR table.
        n_rows: Read only the first ``n_rows`` rows.

    Returns:
        Canonical clonotype frame.

    Raises:
        ValueError: If the signature columns are absent.
    """
    raw = _read_tsv(path, n_rows=n_rows)
    lo = _lower_map(raw.columns)
    count_c = _pick(lo, "read.count")
    nt_c = _pick(lo, "cdr3.nucleotide.sequence")
    aa_c = _pick(lo, "cdr3.amino.acid.sequence")
    v_c = _pick(lo, "v.gene")
    j_c = _pick(lo, "j.gene")
    d_c = _pick(lo, "d.gene")
    if not (count_c and nt_c and aa_c and v_c and j_c):
        raise ValueError(f"not a MiTCR/tcR table; have {raw.columns}")
    return _finalize(raw.select(
        _vdj_expr(v_c).alias(V_CALL), _vdj_expr(d_c).alias(D_CALL),
        _vdj_expr(j_c).alias(J_CALL),
        _upper_expr(nt_c).alias(JUNCTION_NT),
        _unified_expr(aa_c).alias(JUNCTION_AA),
        pl.col(count_c).cast(pl.Float64, strict=False).cast(pl.Int64).alias(COUNT),
    ))


def read_rtcr(path: str | os.PathLike, n_rows: int | None = None) -> pl.DataFrame:
    """Read an RTCR clonotype table (junction aa is re-translated from the nt junction)."""
    raw = _read_tsv(path, n_rows=n_rows)
    lo = _lower_map(raw.columns)
    count_c = _pick(lo, "number of reads")
    v_c = _pick(lo, "v gene")
    j_c = _pick(lo, "j gene")
    nt_c = _pick(lo, "junction nucleotide sequence")
    if not (count_c and v_c and j_c and nt_c):
        raise ValueError(f"not an RTCR table; have {raw.columns}")
    got = raw.select(
        _vdj_expr(v_c).alias(V_CALL), pl.lit(None, pl.Utf8).alias(D_CALL),
        _vdj_expr(j_c).alias(J_CALL), _upper_expr(nt_c).alias(JUNCTION_NT),
        pl.col(count_c).cast(pl.Float64, strict=False).cast(pl.Int64).alias(COUNT),
    )
    return _finalize(got.with_columns(
        _translated(got[JUNCTION_NT]).alias(JUNCTION_AA)))




def _imgt_gene_expr(col: str | None) -> pl.Expr:
    """IMGT/HighV-QUEST ``"Homsap IGHV2-26*01 F"`` -> ``"IGHV2-26"``, over a whole column.

    :func:`_vdj_expr` takes the first tie and strips the allele, leaving ``"Homsap IGHV2-26"``;
    the extract then drops the species prefix and the functionality flag around the gene token.
    """
    if col is None:
        return pl.lit(None, pl.Utf8)
    return _vdj_expr(col).str.extract(r"((?:IG|TR)[A-Z0-9-]+)", 1)


def read_imgt(path: str | os.PathLike, n_rows: int | None = None) -> pl.DataFrame:
    """Read an IMGT/HighV-QUEST ``1_Summary`` table (per-read → collapsed clonotypes)."""
    raw = _read_tsv(path, n_rows=n_rows)
    lo = _lower_map(raw.columns)
    v_c = _pick(lo, "v-gene and allele")
    j_c = _pick(lo, "j-gene and allele")
    d_c = _pick(lo, "d-gene and allele")
    junc_c = _pick(lo, "junction")
    if not (v_c and j_c and junc_c):
        raise ValueError(f"not an IMGT/HighV-QUEST table; have {raw.columns[:8]}…")
    got = (raw.select(
        _imgt_gene_expr(v_c).alias(V_CALL), _imgt_gene_expr(d_c).alias(D_CALL),
        _imgt_gene_expr(j_c).alias(J_CALL), _upper_expr(junc_c).alias(JUNCTION_NT),
        pl.lit(1, pl.Int64).alias(COUNT),   # per-read; _finalize collapses and sums
    ).filter(pl.col(JUNCTION_NT).str.contains(r"^[ATGC]+$")))   # reject empty / N-containing
    return _finalize(got.with_columns(
        _translated(got[JUNCTION_NT]).alias(JUNCTION_AA)))


def read_immunoseq(path: str | os.PathLike, n_rows: int | None = None) -> pl.DataFrame:
    """Read an Adaptive immunoSEQ export (v1 or v2 header dialect, auto-detected).

    The Adaptive nomenclature (``TCRBV29-01``) is converted to IMGT (``TRBV29-1``) with the
    gene→family→family-ties fallback; the CDR3/junction nt is sliced out of the full
    ``rearrangement`` / ``nucleotide`` read via the vIndex + cdr3Length coordinates.
    """
    raw = _read_tsv(path, n_rows=n_rows)
    lo = _lower_map(raw.columns)
    v2 = _pick(lo, "count (templates/reads)") is not None or _pick(lo, "aminoacid") is not None
    if v2:
        count_c = _pick(lo, "count (templates/reads)", "count")
        count2_c = None
        full_c, aa_c, frame_c = _pick(lo, "nucleotide"), _pick(lo, "aminoacid"), _pick(lo, "sequencestatus")
        len_c, idx_c = _pick(lo, "cdr3length"), _pick(lo, "vindex")
        vg, vf, vt = _pick(lo, "vgenename"), _pick(lo, "vfamilyname"), _pick(lo, "vfamilyties")
        dg, df_, dt = _pick(lo, "dgenename"), _pick(lo, "dfamilyname"), _pick(lo, "dfamilyties")
        jg, jf, jt = _pick(lo, "jgenename"), _pick(lo, "jfamilyname"), _pick(lo, "jfamilyties")
    else:
        count_c, count2_c = _pick(lo, "templates"), _pick(lo, "reads")  # templates often "null" → reads
        full_c, aa_c, frame_c = _pick(lo, "rearrangement"), _pick(lo, "amino_acid"), _pick(lo, "frame_type")
        len_c, idx_c = _pick(lo, "cdr3_length"), _pick(lo, "v_index")
        vg, vf, vt = _pick(lo, "v_gene"), _pick(lo, "v_family"), _pick(lo, "v_family_ties")
        dg, df_, dt = _pick(lo, "d_gene"), _pick(lo, "d_family"), _pick(lo, "d_family_ties")
        jg, jf, jt = _pick(lo, "j_gene"), _pick(lo, "j_family"), _pick(lo, "j_family_ties")
    if not (count_c and full_c and len_c and idx_c and vg and jg):
        raise ValueError(f"not an immunoSEQ table; have {raw.columns[:6]}…")

    # The junction is sliced out of the full rearrangement by the vIndex + cdr3Length coordinates,
    # per row -- which `str.slice` takes as expressions, so it stays one pass.
    start = pl.col(idx_c).cast(pl.Int64, strict=False)
    ln = pl.col(len_c).cast(pl.Int64, strict=False)
    nt = (pl.when(start.is_not_null() & ln.is_not_null() & (start >= 0) & (ln > 0))
          .then(pl.col(full_c).cast(pl.Utf8).str.slice(start, ln).str.to_uppercase())
          .otherwise(None))
    count = pl.col(count_c).cast(pl.Float64, strict=False)
    if count2_c:
        # `_to_int` skips non-numeric AND non-positive cells: v1 writes "null" *or* "0" into
        # templates and the real count is in reads. Falling back only on null loses the "0" case.
        alt = pl.col(count2_c).cast(pl.Float64, strict=False)
        count = pl.when(count.is_not_null() & (count > 0)).then(count).otherwise(alt)

    got = raw.select(
        pl.when(nt.str.len_bytes() > 0).then(nt).otherwise(None).alias(JUNCTION_NT),
        (pl.col(frame_c).cast(pl.Utf8).str.strip_chars().str.to_lowercase()
         if frame_c else pl.lit(None, pl.Utf8)).alias("_status"),
        (_unified_expr(aa_c) if aa_c else pl.lit(None, pl.Utf8)).alias("_aa_src"),
        count.cast(pl.Int64).alias(COUNT),
    )
    # Adaptive's own amino acid when it says the read is in-frame, else translate the junction.
    return _finalize(got.with_columns(
        _adaptive_col(raw, vg, vf, vt).alias(V_CALL),
        _adaptive_col(raw, dg, df_, dt).alias(D_CALL),
        _adaptive_col(raw, jg, jf, jt).alias(J_CALL),
        _translated(got[JUNCTION_NT]).alias("_aa_tr"),
    ).with_columns(
        pl.when((pl.col("_status") == "in") & pl.col("_aa_src").is_not_null()
                & (pl.col("_aa_src").str.len_bytes() > 0))
        .then(pl.col("_aa_src")).otherwise(pl.col("_aa_tr")).alias(JUNCTION_AA),
    ).drop("_status", "_aa_src", "_aa_tr"))


def read_vidjil(path: str | os.PathLike, sample_id: int = 0) -> pl.DataFrame:
    """Read a Vidjil ``.vidjil`` JSON file.

    Uses the anchor-inclusive ``seg.junction`` (never ``seg.cdr3``, which excludes the
    anchors); the junction nt is sliced from the clone's full ``sequence`` by the 1-based
    ``junction.start``/``junction.stop`` interval. ``sample_id`` selects the ``reads`` count
    for multi-sample files.
    """
    doc = json.loads(_read_text(path))
    rows = []
    for clone in doc.get("clones", []):
        seg = clone.get("seg")
        if not seg:
            continue
        junction = seg.get("junction")
        if not junction:
            continue
        sequence = clone.get("sequence") or ""
        start, stop = junction.get("start"), junction.get("stop")
        nt = None
        if sequence and start is not None and stop is not None:
            nt = sequence[start - 1:stop].upper() or None  # 1-based inclusive → py slice
        reads = clone.get("reads") or [0]
        # Raise rather than silently fall back to reads[0]: a `sample_id` past the end of a clone's
        # reads list is a user error (wrong index, or 1-based when the format is 0-based), and
        # returning sample 0's counts under a different sample's name is a plausible-looking wrong
        # answer with no error.
        if sample_id >= len(reads):
            raise IndexError(
                f"vidjil sample_id={sample_id} out of range: a clone has {len(reads)} sample(s) "
                f"(valid 0..{len(reads) - 1})"
            )
        cnt = reads[sample_id]
        rows.append({
            V_CALL: extract_vdj(seg.get("5")), D_CALL: extract_vdj(seg.get("4")),
            J_CALL: extract_vdj(seg.get("3")),
            JUNCTION_NT: nt, JUNCTION_AA: to_unified_cdr3aa(junction.get("aa")),
            COUNT: _to_int(cnt),
        })
    # The loop stays: this is a JSON document, not a table, so there is no column to express the
    # work over -- the per-clone Python is the parsing itself. It hands `_finalize` a frame with a
    # declared schema so an all-null D call cannot come out as dtype Null.
    return _finalize(pl.DataFrame(rows, schema={
        V_CALL: pl.Utf8, D_CALL: pl.Utf8, J_CALL: pl.Utf8,
        JUNCTION_NT: pl.Utf8, JUNCTION_AA: pl.Utf8, COUNT: pl.Int64,
    }))


def read_trust4(path: str | os.PathLike, n_rows: int | None = None) -> pl.DataFrame:
    """Read a TRUST4 clonotype report (``*_report.tsv``).

    The TRUST4 report header is
    ``#count  frequency  CDR3nt  CDR3aa  V  D  J  C  cid  cid_full_length``. TRUST4's CDR3
    spans the conserved Cys104 … Phe/Trp118 anchors (anchors **included**), i.e. it is the
    AIRR junction, so ``CDR3nt``/``CDR3aa`` map straight to ``junction_nt``/``junction_aa``
    (the ``_`` stop / ``?`` ambiguous-N markers are collapsed by :func:`to_unified_cdr3aa`).
    ``V``/``D``/``J``/``C`` keep the first allele's gene (``*`` → missing); the C column is
    the BCR isotype (or the TCR constant gene) when the constant region was captured. Rows
    whose CDR3 nt is not clean ``ACGT`` (TRUST4 ``partial`` / ``out_of_frame`` / N-containing)
    are dropped.
    """
    raw = _read_tsv(path, n_rows=n_rows)
    lo = {c.lower().lstrip("#"): c for c in raw.columns}  # first column is ``#count``
    count_c = _pick(lo, "count")
    nt_c, aa_c = _pick(lo, "cdr3nt"), _pick(lo, "cdr3aa")
    v_c, d_c, j_c, c_c = _pick(lo, "v"), _pick(lo, "d"), _pick(lo, "j"), _pick(lo, "c")
    if not (count_c and nt_c and aa_c and v_c and j_c):
        raise ValueError(f"not a TRUST4 report (need count / CDR3nt+aa / V,J); have {raw.columns}")
    return _finalize(raw.select(
        _vdj_expr(v_c).alias(V_CALL), _vdj_expr(d_c).alias(D_CALL),
        _vdj_expr(j_c).alias(J_CALL), _vdj_expr(c_c).alias(C_CALL),
        _upper_expr(nt_c).alias(JUNCTION_NT), _unified_expr(aa_c).alias(JUNCTION_AA),
        pl.col(count_c).cast(pl.Float64, strict=False).cast(pl.Int64).alias(COUNT),
    # skip TRUST4 partial / out-of-frame / N-containing CDR3s
    ).filter(pl.col(JUNCTION_NT).str.contains(r"^[ATGC]+$")))


def read_arda(path: str | os.PathLike, n_rows: int | None = None) -> pl.DataFrame:
    """Read arda's AIRR annotation output (per-sequence ``*.airr.tsv`` or ``clones.tsv``).

    arda (AIRR annotation + markup repair) writes standard AIRR Rearrangement column names,
    so this delegates to :func:`~vdjtools.io.read.read_airr` — which maps
    ``v_call``/``d_call``/``j_call``/``c_call`` and ``junction``/``junction_aa`` and collapses
    reads to clonotypes — then nulls the literal ``""`` arda emits for an empty gene call.
    arda's extra columns (``d2_call``, ``c_class``, ``mmseqs2_*``) are ignored.
    """
    df = read_airr(path, n_rows=n_rows)
    return df.with_columns(
        pl.when(pl.col(c).str.strip_chars('"') == "").then(None).otherwise(pl.col(c)).alias(c)
        for c in (D_CALL, C_CALL) if c in df.columns
    )
