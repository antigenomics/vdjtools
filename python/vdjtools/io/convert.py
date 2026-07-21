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

import gzip
import json
import os
import re
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


def _adaptive_to_imgt(field: str | None) -> str | None:
    """Adaptive immunoSEQ → IMGT: ``extract_vdj`` then ``TCR→TR`` and drop zero-padding.

    ``"TCRBV29-01" → "TRBV29-1"``; ``"unresolved" → None``. Port of
    ``CommonUtil.extractVDJImmunoSeq`` (single-value form).
    """
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


def _finalize(rows: list[dict]) -> pl.DataFrame:
    """Filter bad rows, collapse to unique clonotypes (summed counts), coerce to canonical."""
    if not rows:
        return schema.add_locus(schema.normalize(pl.DataFrame(schema={c: pl.Utf8 for c in
                                (V_CALL, D_CALL, J_CALL, JUNCTION_AA, JUNCTION_NT)})))
    df = pl.DataFrame(rows)
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
    rows = []
    for r in raw.iter_rows(named=True):
        cnt = r[count_c]
        rows.append({
            V_CALL: extract_vdj(r[v_c]), D_CALL: extract_vdj(r[d_c]) if d_c else None,
            J_CALL: extract_vdj(r[j_c]), C_CALL: extract_vdj(r[c_c]) if c_c else None,
            JUNCTION_NT: (r[nt_c] or "").upper() or None,
            JUNCTION_AA: r[aa_c],  # MiXcr aa is milib-based — kept verbatim (no unify)
            COUNT: _to_int(cnt),
        })
    return _finalize(rows)


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
    rows = [{
        V_CALL: extract_vdj(r[v_c]), D_CALL: extract_vdj(r[d_c]) if d_c else None,
        J_CALL: extract_vdj(r[j_c]),
        JUNCTION_NT: (r[nt_c] or "").upper() or None,
        JUNCTION_AA: to_unified_cdr3aa(r[aa_c]),
        COUNT: _to_int(r[count_c]),
    } for r in raw.iter_rows(named=True)]
    return _finalize(rows)


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
    rows = [{
        V_CALL: extract_vdj(r[v_c]), D_CALL: extract_vdj(r[d_c]) if d_c else None,
        J_CALL: extract_vdj(r[j_c]),
        JUNCTION_NT: (r[nt_c] or "").upper() or None,
        JUNCTION_AA: to_unified_cdr3aa(r[aa_c]),
        COUNT: _to_int(r[count_c]),
    } for r in raw.iter_rows(named=True)]
    return _finalize(rows)


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
    rows = []
    for r in raw.iter_rows(named=True):
        nt = (r[nt_c] or "").upper() or None
        rows.append({
            V_CALL: extract_vdj(r[v_c]), D_CALL: None, J_CALL: extract_vdj(r[j_c]),
            JUNCTION_NT: nt,
            JUNCTION_AA: to_unified_cdr3aa(translate(nt)) if nt else None,
            COUNT: _to_int(r[count_c]),
        })
    return _finalize(rows)


_IMGT_GENE = re.compile(r"(?:IG|TR)[A-Z0-9-]+")
_ATGC_ONLY = re.compile(r"^[ATGC]+$")


def _imgt_gene(field: str | None) -> str | None:
    """IMGT/HighV-QUEST ``"Homsap IGHV2-26*01 F"`` → ``"IGHV2-26"`` (strip species/allele/flag)."""
    gene = extract_vdj(field)  # first tie, allele stripped -> "Homsap IGHV2-26"
    if gene is None:
        return None
    m = _IMGT_GENE.search(gene)
    return m.group(0) if m else None


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
    rows = []
    for r in raw.iter_rows(named=True):
        nt = (r[junc_c] or "").upper()
        if not _ATGC_ONLY.match(nt):  # reject empty / N-containing junctions
            continue
        rows.append({
            V_CALL: _imgt_gene(r[v_c]), D_CALL: _imgt_gene(r[d_c]) if d_c else None,
            J_CALL: _imgt_gene(r[j_c]),
            JUNCTION_NT: nt, JUNCTION_AA: to_unified_cdr3aa(translate(nt)),
            COUNT: 1,  # per-read output; _finalize collapses identical junctions and sums
        })
    return _finalize(rows)


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

    def _cell(r, c):
        return r[c] if c else None

    rows = []
    for r in raw.iter_rows(named=True):
        full, vidx, clen = r[full_c], r[idx_c], r[len_c]
        nt = None
        if full and vidx not in (None, "") and clen not in (None, ""):
            start, ln = int(vidx), int(clen)
            if start >= 0 and ln > 0:
                nt = full[start:start + ln].upper() or None
        status = (_cell(r, frame_c) or "").strip().lower()
        aa_src = _cell(r, aa_c)
        if status == "in" and aa_src:
            junc_aa = to_unified_cdr3aa(aa_src)
        else:
            junc_aa = to_unified_cdr3aa(translate(nt)) if nt else None
        rows.append({
            V_CALL: _adaptive_call(_cell(r, vg), _cell(r, vf), _cell(r, vt)),
            D_CALL: _adaptive_call(_cell(r, dg), _cell(r, df_), _cell(r, dt)),
            J_CALL: _adaptive_call(_cell(r, jg), _cell(r, jf), _cell(r, jt)),
            JUNCTION_NT: nt, JUNCTION_AA: junc_aa,
            COUNT: _to_int(r[count_c], r[count2_c] if count2_c else None),
        })
    return _finalize(rows)


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
    return _finalize(rows)


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
    rows = []
    for r in raw.iter_rows(named=True):
        nt = (r[nt_c] or "").upper()
        if not _ATGC_ONLY.match(nt):  # skip TRUST4 partial / out-of-frame / N-containing CDR3s
            continue
        rows.append({
            V_CALL: extract_vdj(r[v_c]), D_CALL: extract_vdj(r[d_c]) if d_c else None,
            J_CALL: extract_vdj(r[j_c]), C_CALL: extract_vdj(r[c_c]) if c_c else None,
            JUNCTION_NT: nt, JUNCTION_AA: to_unified_cdr3aa(r[aa_c]),
            COUNT: _to_int(r[count_c]),
        })
    return _finalize(rows)


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
