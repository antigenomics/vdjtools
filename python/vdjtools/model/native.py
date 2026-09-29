"""Bridge from the polars :class:`Model` to the native ``_core`` hot loops.

``pack`` reconstructs the model's dense arrays (gene choice, deletion, insertion, dinucleotide)
and germline cut segments into a C++ :class:`PackedModel`; ``pgen_nt`` calls the native Pgen.
The result matches the pure-Python reference (and OLGA) exactly — the native path is just faster.
"""
from __future__ import annotations

from .model import Model
from .pgen import prepare

_NT2NUM = {"A": 0, "C": 1, "G": 2, "T": 3}
_pack_cache: dict[int, tuple] = {}


def _encode(seq: str) -> list[int]:
    return [_NT2NUM[c] for c in seq]


def _gene_idx(idx_of: dict[str, int], name: str | None, kind: str) -> int:
    """Resolve a V/J **allele** name to its model index; ``None`` marginalizes (index ``-1``).

    Raises on an unknown name rather than falling back to ``-1``: ``-1`` means "marginalize over
    every allele", so a silent fallback turns a mis-typed or gene-level call into a *different,
    larger* Pgen with no error — e.g. ``"TRBV9"`` (gene) returned the V/J-agnostic value, 2.38x
    the true ``"TRBV9*01"`` Pgen.
    """
    if not name:
        return -1
    idx = idx_of.get(name)
    if idx is not None:
        return idx
    alleles = sorted(a for a in idx_of if a.split("*")[0] == name)
    if alleles:
        raise KeyError(
            f"{kind} {name!r} is a gene name; the model is keyed by allele. Pass one of "
            f"{alleles}, or None to marginalize over all {kind}."
        )
    raise KeyError(
        f"{kind} {name!r} is not in the model ({len(idx_of)} alleles, e.g. "
        f"{sorted(idx_of)[0]!r}). Pass a known allele, or None to marginalize over all {kind}."
    )


def gene_to_allele(model: Model) -> dict[str, str]:
    """Map each V/J **gene** the model carries to a representative allele of that gene.

    Real repertoires and every aligner report genes (``TRBV10-3``) while a model is keyed by
    allele (``TRBV10-3*01``), and :func:`pgen_aa` deliberately raises on a gene name rather than
    silently marginalising over every allele — that fallback once returned a Pgen 2.38x too high
    with no error. So the gene has to be resolved to a concrete allele deliberately and visibly.

    The representative is the lowest-numbered allele present, i.e. ``*01`` wherever the model has
    it. Alleles of one gene share the CDR3-region germline in all but rare cases, so this is the
    conventional reading of a gene-level call — but it IS a choice, which is why the functions
    using it expose ``resolve_genes=False`` to refuse it instead.

    Args:
        model: A recombination :class:`Model`.

    Returns:
        ``{gene: allele}`` over both V and J (their gene names are disjoint).
    """
    _pm, vi, ji = pack(model)
    out: dict[str, str] = {}
    for idx_of in (vi, ji):
        for allele in idx_of:
            gene = allele.split("*")[0]
            if gene not in out or allele < out[gene]:
                out[gene] = allele
    return out


def _resolve_vj(model: Model, v: str | None, j: str | None,
                resolve_genes: bool) -> tuple[str | None, str | None]:
    """``(v, j)`` with a gene-level call replaced by its representative allele, when asked."""
    if not resolve_genes:
        return v, j
    alias = gene_to_allele(model)
    return (alias.get(v, v) if v else v), (alias.get(j, j) if j else j)


def _del_dense(pdel: dict, idx_of: dict, maxpal: int) -> tuple[list[float], int]:
    """{(allele, ndel): p} -> (flat [n_allele * nbins], nbins), index = ndel + maxpal."""
    nbins = max((n + maxpal for (_a, n) in pdel), default=0) + 1
    arr = [0.0] * (len(idx_of) * nbins)
    for (a, n), p in pdel.items():
        arr[idx_of[a] * nbins + (n + maxpal)] = p
    return arr, nbins


def _del_dense_d(pdel: dict, idx_of: dict, maxdl: int, maxdr: int) -> tuple[list[float], int, int]:
    n5 = max((a + maxdl for (_d, a, _b) in pdel), default=0) + 1
    n3 = max((b + maxdr for (_d, _a, b) in pdel), default=0) + 1
    arr = [0.0] * (len(idx_of) * n5 * n3)
    for (d, a, b), p in pdel.items():
        arr[(idx_of[d] * n5 + (a + maxdl)) * n3 + (b + maxdr)] = p
    return arr, n5, n3


def _del_dense_d_fixed(pdel: dict, idx_of: dict, maxdl: int, maxdr: int, n5: int, n3: int) -> list[float]:
    """Pack a D-deletion dict into an explicit ``[nD*n5*n3]`` grid (for the second D, which must
    share the first D's ``nbins`` so the C++ index arithmetic matches)."""
    arr = [0.0] * (len(idx_of) * n5 * n3)
    for (d, a, b), p in pdel.items():
        i5, i3 = a + maxdl, b + maxdr
        if 0 <= i5 < n5 and 0 <= i3 < n3:
            arr[(idx_of[d] * n5 + i5) * n3 + i3] = p
    return arr


def pack(model: Model):
    """Build (and cache) the native :class:`PackedModel` for ``model``.

    The native nt Pgen, aa Pgen (incl. Hamming-1 and v/j-agnostic), and the EM E-step
    (:func:`~vdjtools.model.infer.infer_native`) all support tandem-D (``n_D=2``).

    The cache is keyed by ``id(model)`` but stores the model reference and verifies identity on hit:
    CPython recycles object ids after GC, so a bare-id cache could return a stale :class:`PackedModel`
    for a *different* model that reused a freed id (e.g. running TRB then TRD EM in one process — the
    stale TRB pack has a different gene count and crashes the M-step). Keeping the ref also pins the id.
    """
    key = id(model)
    hit = _pack_cache.get(key)
    if hit is not None and hit[3] is model:
        return hit[:3]

    from .._core import PackedModel

    prep = prepare(model)
    vdj = model.chain_type == "VDJ"
    v_alleles = model.genomic["genes_v"]["v_allele"].to_list()
    j_alleles = model.genomic["genes_j"]["j_allele"].to_list()
    d_alleles = model.genomic["genes_d"]["d_allele"].to_list() if vdj else []
    vi = {a: i for i, a in enumerate(v_alleles)}
    ji = {a: i for i, a in enumerate(j_alleles)}
    di = {a: i for i, a in enumerate(d_alleles)}

    pm = PackedModel()
    pm.vdj = vdj
    pm.maxpal_v3 = prep.maxpal["v_3"]
    pm.maxpal_j5 = prep.maxpal["j_5"]
    pm.cut_v = [_encode(prep.cut["v"][a]) for a in v_alleles]
    pm.cut_j = [_encode(prep.cut["j"][a]) for a in j_alleles]
    pm.func_v = [vi[a] for a in prep.functional_v]
    pm.func_j = [ji[a] for a in prep.functional_j]
    pm.pv = [float(prep.p_v.get(a, 0.0)) for a in v_alleles]
    pm.del_v, pm.nbins_v = _del_dense(prep.p_del["v"], vi, prep.maxpal["v_3"])
    pm.del_j, pm.nbins_j = _del_dense(prep.p_del["j"], ji, prep.maxpal["j_5"])

    if vdj:
        pm.maxpal_d5 = prep.maxpal["d_5"]
        pm.maxpal_d3 = prep.maxpal["d_3"]
        pm.cut_d = [_encode(prep.cut["d"][a]) for a in d_alleles]
        pm.func_d = [di[a] for a in prep.functional_d]
        pm.pj = [float(prep.p_j.get(a, 0.0)) for a in j_alleles]
        pm.pd_given_j = [float(prep.p_d_given_j.get((j, d), 0.0)) for j in j_alleles for d in d_alleles]
        pm.del_d, pm.nbins_d5, pm.nbins_d3 = _del_dense_d(prep.p_del["d"], di, prep.maxpal["d_5"], prep.maxpal["d_3"])
        pm.ins_vd = prep.p_ins["vd"].tolist()
        pm.ins_dj = prep.p_ins["dj"].tolist()
        pm.R_vd = prep.R["vd"].reshape(-1).tolist()
        pm.R_dj = prep.R["dj"].reshape(-1).tolist()
        pm.bias_vd = prep.bias["vd"].tolist()
        pm.bias_dj = prep.bias["dj"].tolist()
        # D-D (tandem) extension — populated only when the model declares it.
        pm.p_nd1 = float(prep.p_nd.get(0, 0.0) + prep.p_nd.get(1, 0.0))
        pm.p_nd2 = float(prep.p_nd.get(2, 0.0))
        pm.dd = bool(prep.p_d2_given_d1)
        if pm.dd:
            pm.pd2_given_d1 = [float(prep.p_d2_given_d1.get((d1, d2), 0.0)) for d1 in d_alleles for d2 in d_alleles]
            pm.del_d2 = _del_dense_d_fixed(prep.p_del_d2, di, prep.maxpal["d_5"], prep.maxpal["d_3"], pm.nbins_d5, pm.nbins_d3)
            pm.ins_dd = prep.p_ins["dd"].tolist()
            pm.R_dd = prep.R["dd"].reshape(-1).tolist()
            pm.bias_dd = prep.bias["dd"].tolist()
    else:
        pm.pjv = [float(prep.p_j.get((v, j), 0.0)) for v in v_alleles for j in j_alleles]
        pm.ins_vj = prep.p_ins["vj"].tolist()
        pm.R_vj = prep.R["vj"].reshape(-1).tolist()
        pm.bias_vj = prep.bias["vj"].tolist()

    _pack_cache[key] = (pm, vi, ji, model)
    return pm, vi, ji


def pgen_nt(model: Model, cdr3_nt: str, v: str | None = None, j: str | None = None) -> float:
    """Native nucleotide Pgen — same result as :func:`vdjtools.model.pgen.pgen_nt`, faster."""
    from .._core import pgen_nt as _pgen_nt

    pm, vi, ji = pack(model)
    return _pgen_nt(pm, _encode(cdr3_nt.upper()),
                    _gene_idx(vi, v, "V"), _gene_idx(ji, j, "J"))


def pgen_aa(
    model: Model,
    cdr3_aa: str,
    v: str | None = None,
    j: str | None = None,
    mismatches: int = 0,
) -> float:
    """Native amino-acid Pgen — same result as :func:`vdjtools.model.pgen.pgen_aa`, much faster.

    Args:
        model: A recombination :class:`Model`.
        cdr3_aa: The junction/CDR3 amino-acid sequence (Cys → Phe/Trp inclusive).
        v: V **allele** to condition on (e.g. ``"TRBV9*01"``), or ``None`` to marginalize over all
            V (V-agnostic). A gene-level name (``"TRBV9"``) or any unknown allele raises
            :class:`KeyError` — it must not silently degrade to the V-agnostic value.
        j: J allele to condition on (e.g. ``"TRBJ2-3*01"``), or ``None`` to marginalize (as ``v``).
        mismatches: ``0`` for the exact sequence; ``1`` to also sum the Pgen of every
            amino-acid sequence within Hamming distance 1 (one substitution) — the total
            probability mass in the 1-mismatch ball, computed natively far faster than OLGA's
            per-neighbour approach.

    Returns:
        Pgen as a float.
    """
    from .._core import pgen_aa as _pgen_aa
    from .._core import pgen_aa_hamming1 as _pgen_aa_h1

    pm, vi, ji = pack(model)
    vidx = _gene_idx(vi, v, "V")
    jidx = _gene_idx(ji, j, "J")
    if mismatches == 0:
        return _pgen_aa(pm, cdr3_aa.upper(), vidx, jidx)
    if mismatches == 1:
        return _pgen_aa_h1(pm, cdr3_aa.upper(), vidx, jidx)
    raise ValueError("mismatches must be 0 or 1")


def pgen_aa_degenerate(
    model: Model,
    allowed: list[str],
    v: str | None = None,
    j: str | None = None,
) -> float:
    """Total Pgen of every CDR3 matching a motif — a per-position set of permitted residues.

    The same masked transfer-matrix DP :func:`pgen_aa` runs (that is the all-singletons case) and
    :func:`pgen_aa` with ``mismatches=1`` runs per position (the one-wildcard case), with the
    residue set at each position free. Cost is one DP pass whatever the sets contain: the sequences
    the motif matches are summed over, never enumerated. This is the exact Pgen of a V/J/length-
    pinned motif such as a VDJdb cluster PWM thresholded per position.

    Args:
        model: A recombination :class:`Model`.
        allowed: One item per amino-acid position of the junction/CDR3 (so ``len(allowed)`` is the
            junction length). Each item is a string of the residues permitted there: ``"C"`` pins
            one, ``"ILVF"`` allows a subset, and ``""`` or a string containing ``"X"`` is a
            **wildcard** (any of the 20 amino acids). ``list(seq)`` reproduces
            :func:`pgen_aa` exactly. A character the genetic code does not name raises
            :class:`ValueError` — an empty codon mask would score a silent ``0.0``.
        v: V **allele** to condition on, or ``None`` to marginalize — as :func:`pgen_aa`, including
            the :class:`KeyError` on a gene-level or unknown name.
        j: J allele to condition on, or ``None`` to marginalize (as ``v``).

    Returns:
        Pgen as a float — the summed Pgen of all matching sequences.

    Note:
        ``X`` means wildcard **here only**. :func:`pgen_aa` matches residues by exact character
        against the genetic code, so ``pgen_aa(m, "CASSLAPGATNEKLXF")`` is ``0.0``, not a
        wildcard query. Pass ``list("CASSLAPGATNEKLXF")`` to this function to get the latter.
    """
    from .._core import pgen_aa_degenerate as _deg

    pm, vi, ji = pack(model)
    return _deg(pm, [str(a).upper() for a in allowed],
                _gene_idx(vi, v, "V"), _gene_idx(ji, j, "J"))


def pgen_aa_degenerate_batch(
    model: Model,
    allowed: list[list[str]],
    v: list[str | None] | None = None,
    j: list[str | None] | None = None,
    threads: int = 0,
) -> list[float]:
    """Batch :func:`pgen_aa_degenerate` over many motifs, parallelized across queries.

    Mirrors :func:`pgen_aa_batch`: the GIL is released and the queries are partitioned across
    worker threads, and the result is bitwise-identical to the serial per-motif computation for
    any ``threads``.

    Args:
        model: A recombination :class:`Model`.
        allowed: One per-position residue-set list per motif (see :func:`pgen_aa_degenerate`).
        v: Optional per-motif V **alleles** (same length as ``allowed``); ``None`` marginalizes
            over all V for every motif. Individual entries may be ``None``.
        j: Optional per-motif J alleles (as ``v``).
        threads: Worker threads; ``0`` = auto (``hardware_concurrency - 2``). Batches under 64
            motifs run single-threaded.

    Returns:
        Per-motif Pgen in input order.
    """
    from .._core import pgen_aa_degenerate_batch as _batch

    pm, vi, ji = pack(model)
    sets = [[str(a).upper() for a in one] for one in allowed]
    v_idxs = [_gene_idx(vi, x, "V") for x in v] if v is not None else []
    j_idxs = [_gene_idx(ji, x, "J") for x in j] if j is not None else []
    if v_idxs and len(v_idxs) != len(sets):
        raise ValueError("v must have the same length as allowed")
    if j_idxs and len(j_idxs) != len(sets):
        raise ValueError("j must have the same length as allowed")
    return _batch(pm, sets, v_idxs, j_idxs, threads)


def pgen_aa_batch(
    model: Model,
    cdr3_aas: list[str],
    v: list[str | None] | None = None,
    j: list[str | None] | None = None,
    mismatches: int = 0,
    threads: int = 0,
) -> list[float]:
    """Batch amino-acid Pgen over many CDR3s, parallelized across sequences in native code.

    Computes the same value as calling :func:`pgen_aa` per sequence, but releases the GIL and
    partitions the sequences across worker threads — the clean, exact speedup for the real
    workload (Pgen / 1-mismatch matching over many clonotypes). The result is bitwise-identical
    to the serial per-sequence computation for any ``threads``.

    Args:
        model: A recombination :class:`Model`.
        cdr3_aas: Junction/CDR3 amino-acid sequences.
        v: Optional per-sequence V **alleles** to condition on (same length as ``cdr3_aas``);
            ``None`` marginalises over all V for every sequence. Individual entries may be
            ``None``. An unknown or gene-level name raises :class:`KeyError` (see :func:`pgen_aa`).
        j: Optional per-sequence J alleles (as ``v``).
        mismatches: ``0`` for exact Pgen, ``1`` for the Hamming-1 ball (as :func:`pgen_aa`).
        threads: Worker threads; ``0`` = auto (``hardware_concurrency - 2``). Batches under 64
            sequences run single-threaded.

    Returns:
        Per-sequence Pgen in input order.
    """
    if mismatches not in (0, 1):
        raise ValueError("mismatches must be 0 or 1")
    from .._core import pgen_aa_batch as _batch

    pm, vi, ji = pack(model)
    seqs = [s.upper() for s in cdr3_aas]
    v_idxs = [_gene_idx(vi, x, "V") for x in v] if v is not None else []
    j_idxs = [_gene_idx(ji, x, "J") for x in j] if j is not None else []
    if v_idxs and len(v_idxs) != len(seqs):
        raise ValueError("v must have the same length as cdr3_aas")
    if j_idxs and len(j_idxs) != len(seqs):
        raise ValueError("j must have the same length as cdr3_aas")
    return _batch(pm, seqs, v_idxs, j_idxs, mismatches, threads)


def pgen_nt_batch(model: Model, cdr3_nts: list[str], v: list | None = None,
                  j: list | None = None, threads: int = 0) -> list[float]:
    """Nucleotide Pgen over many sequences, parallelized natively across the batch.

    The counterpart of :func:`pgen_aa_batch`, and the entry point to use on more than one
    sequence. Calling :func:`pgen_nt` in a loop pays the per-call Python each time, and wrapping
    that loop in a thread pool dispatches one task per sequence -- which is what this replaces.

    Args:
        model: A recombination :class:`Model`.
        cdr3_nts: Junction/CDR3 nucleotide sequences.
        v: Optional per-sequence V **alleles** (same length as ``cdr3_nts``); ``None`` marginalises
            over all V for every sequence, and individual entries may be ``None``. A gene-level
            name raises, as in :func:`pgen_nt`.
        j: Optional per-sequence J alleles (as ``v``).
        threads: **Kernel threads**, not worker processes; ``0`` = auto
            (``hardware_concurrency - 2``). Batches under 64 sequences stay single-threaded so the
            result is bitwise-identical to a serial run.

    Returns:
        Per-sequence Pgen in input order.

    Raises:
        KeyError: If any ``v``/``j`` entry names no allele the model carries.
        ValueError: If ``v``/``j`` has the wrong length, or a sequence has a non-ACGT base.
    """
    from .._core import pgen_nt_batch as _batch

    pm, vi, ji = pack(model)
    seqs = [s.upper() for s in cdr3_nts]
    v_idxs = [_gene_idx(vi, x, "V") for x in v] if v is not None else []
    j_idxs = [_gene_idx(ji, x, "J") for x in j] if j is not None else []
    if v_idxs and len(v_idxs) != len(seqs):
        raise ValueError("v must have the same length as cdr3_nts")
    if j_idxs and len(j_idxs) != len(seqs):
        raise ValueError("j must have the same length as cdr3_nts")
    return _batch(pm, seqs, v_idxs, j_idxs, threads)


def best_aa_scenarios(model: Model, cdr3_aa: str, v: str | None = None, j: str | None = None,
                      k: int = 8, *, resolve_genes: bool = True) -> list[tuple]:
    """Top-``k`` recombination scenarios for an amino-acid CDR3, best first.

    The argmax counterpart of :func:`pgen_aa`, over the same Pi_L*Pi_R transfer matrix: ``max`` in
    place of the sums, with the winning ``(V, delV)`` and ``(J, delJ)`` carried through the state.
    Because the DP marginalizes over V and J at no extra cost, leaving both unspecified is barely
    slower than pinning them (0.26 vs 0.23 ms on human TRB).

    Args:
        model: The recombination model.
        cdr3_aa: CDR3 amino-acid sequence (conserved Cys -> conserved Phe/Trp inclusive).
        v, j: Optional V/J names to condition on. ``None`` marginalizes. A gene-level name is
            resolved to a representative allele unless ``resolve_genes=False``; a name the model
            has no gene for raises :class:`KeyError` naming it.
        k: How many scenarios to return.
        resolve_genes: Resolve a **gene**-level call (``TRBV10-3``) to a representative model
            allele (``TRBV10-3*01``) before the DP — see :func:`gene_to_allele`. Default ``True``,
            because real V/J calls are gene-level and without it every one of them raises. Set
            ``False`` to accept exact allele names only.

    Returns:
        ``[(w, v_allele, len_v, j_allele, len_j, d_allele | None, idx5, idx3, pos)]``, descending by
        ``w``. ``len_v``/``len_j`` are the nucleotides each germline contributes to the CDR3, so the
        boundaries are ``v_end = len_v`` and ``j_start = 3 * len(cdr3_aa) - len_j``;
        ``idx5``/``idx3`` are the D's 5'/3' trims and ``pos`` where its contribution starts.

        **Empty** only when the DP can explain nothing: a residue outside the genetic code (which
        has no codons), or no rearrangement reachable under the pinned V/J. It is never empty
        because of a naming mismatch — that raises.

    Raises:
        KeyError: If ``v`` or ``j`` names no gene the model carries (with ``resolve_genes``), or no
            exact allele (without it).

    Note:
        This is the argmax of a *probability* model, not an aligner. Measured on 25,000 real human
        TRB clonotypes against the observed nucleotide markup, the top-1 scenario places ``v_end``
        exactly 59.8% of the time and ``j_start`` 90.3%, against 73.3% / 97.7% for germline
        alignment (``arda.cdr3fix``) at a twelfth of the cost — many junctions have several
        near-equally-probable boundaries, so the modal scenario is often not the one that happened.
        Reach for this when you want what an aligner cannot give: the *alternatives* with their
        probabilities, D geometry, or a V/J named by marginalising over the model. When the V/J call
        is already known and you only want the boundary, align instead.
    """
    from .._core import best_aa_scenarios as _best

    pm, vi, ji = pack(model)
    v, j = _resolve_vj(model, v, j, resolve_genes)
    d_alleles = (model.genomic["genes_d"]["d_allele"].to_list()
                 if model.chain_type == "VDJ" else [])
    iv = {i: a for a, i in vi.items()}
    ij = {i: a for a, i in ji.items()}
    got = _best(pm, cdr3_aa.upper(), _gene_idx(vi, v, "V"), _gene_idx(ji, j, "J"), k)
    return [(s.w, iv[s.v], s.len_v, ij[s.j], s.len_j,
             d_alleles[s.d] if s.d >= 0 else None, s.idx5, s.idx3, s.pos) for s in got]


def best_aa_scenarios_batch(
    model: Model,
    cdr3_aas: list[str],
    v: list[str | None] | None = None,
    j: list[str | None] | None = None,
    k: int = 8,
    threads: int = 0,
    *,
    resolve_genes: bool = True,
):
    """Batch :func:`best_aa_scenarios` over many CDR3s, parallelized across sequences natively.

    Mirrors :func:`pgen_aa_batch`: the GIL is released and the queries are partitioned across
    worker threads, and the result is identical to the serial per-sequence calls for any
    ``threads``. This is the entry point to use on a clonotype table — the per-row Python loop,
    not the DP, is most of the cost of calling :func:`best_aa_scenarios` 25,000 times.

    Args:
        model: A recombination :class:`Model`.
        cdr3_aas: Junction/CDR3 amino-acid sequences.
        v: Optional per-sequence V names (same length as ``cdr3_aas``); ``None`` marginalizes over
            all V for every sequence, and individual entries may be ``None``.
        j: Optional per-sequence J names (as ``v``).
        k: Scenarios per sequence.
        threads: Worker threads; ``0`` = auto (``hardware_concurrency - 2``). Batches under 64
            sequences run single-threaded.
        resolve_genes: As :func:`best_aa_scenarios`.

    Returns:
        A :class:`polars.DataFrame` with one row per **scenario**, ordered by ``row`` then ``rank``:
        ``row`` (index into ``cdr3_aas``), ``rank`` (0-based within that row's top-``k``), ``w``,
        ``v_call``, ``len_v``, ``j_call``, ``len_j``, ``d_call``, ``idx5``, ``idx3``, ``pos``.
        A sequence the DP cannot explain contributes **no rows**, so an absent ``row`` value is how
        a declined query shows up; it is never a silently marginalized scenario.

    Raises:
        KeyError: If any ``v``/``j`` entry names no gene the model carries (see
            :func:`best_aa_scenarios`).
        ValueError: If ``v`` or ``j`` is given with a length other than ``len(cdr3_aas)``.
    """
    import numpy as np
    import polars as pl

    from .._core import best_aa_scenarios_batch as _batch

    pm, vi, ji = pack(model)
    seqs = [s.upper() for s in cdr3_aas]
    if v is not None and len(v) != len(seqs):
        raise ValueError("v must have the same length as cdr3_aas")
    if j is not None and len(j) != len(seqs):
        raise ValueError("j must have the same length as cdr3_aas")
    alias = gene_to_allele(model) if resolve_genes else {}
    v_idxs = ([_gene_idx(vi, alias.get(x, x) if x else x, "V") for x in v]
              if v is not None else [])
    j_idxs = ([_gene_idx(ji, alias.get(x, x) if x else x, "J") for x in j]
              if j is not None else [])
    cols = _batch(pm, seqs, v_idxs, j_idxs, k, threads)

    # Index -> allele name by fancy-indexing, so naming k scenarios per row stays vectorized.
    names_v = np.array([a for a, _i in sorted(vi.items(), key=lambda kv: kv[1])], dtype=object)
    names_j = np.array([a for a, _i in sorted(ji.items(), key=lambda kv: kv[1])], dtype=object)
    names_d = np.array(model.genomic["genes_d"]["d_allele"].to_list()
                       if model.chain_type == "VDJ" else [], dtype=object)
    d = cols["d"]
    return pl.DataFrame({
        "row": cols["row"], "rank": cols["rank"], "w": cols["w"],
        "v_call": names_v[cols["v"]], "len_v": cols["len_v"],
        "j_call": names_j[cols["j"]], "len_j": cols["len_j"],
        "d_call": (np.where(d >= 0, names_d[np.maximum(d, 0)], None)
                   if names_d.size else np.full(d.shape, None, dtype=object)),
        "idx5": cols["idx5"], "idx3": cols["idx3"], "pos": cols["pos"],
    })


def _alts(spec, alias: dict, idx_of: dict, kind: str) -> list[int]:
    """A ``v=``/``j=`` entry -> the gene indices to search.

    Three input modes, because real annotation tables carry all three: ``None`` (nothing known, so
    marginalize over the segment), one name, or **several** -- a list, or the comma-separated string
    an AIRR ``v_call`` holds when the aligner could not choose. The multi case is not a pin: every
    listed allele is searched and the model picks the most plausible.

    An **empty** candidate set is a fourth thing and not the same as ``None``: it names no allele,
    so it explains nothing and the row is declined.
    """
    if spec is None:
        return [-1]
    names = [s.strip() for s in spec.split(",")] if isinstance(spec, str) else list(spec)
    # An EMPTY candidate set is not ``None``: it explains nothing, and the row is declined. Folding
    # the two together would answer a call that named no allele by marginalizing over all of them.
    return [_gene_idx(idx_of, alias.get(n, n), kind) for n in names if n]


def infer_nt_many(model: Model, cdr3_aas: list, v: list | None = None, j: list | None = None,
                  k: int = 8, threads: int = 0, *, resolve_genes: bool = True) -> dict:
    """:func:`vdjtools.model.infer_nt` over many CDR3s, entirely in the native layer.

    The whole pipeline -- scenario search, codon reconstruction, exact marginal re-score -- runs in
    one call with the GIL released, which is what makes ``threads`` worth anything: the
    reconstruction used to be Python and held the GIL for 35% of a human TRB row (#181).

    Args:
        model: A recombination :class:`Model`.
        cdr3_aas: Junction amino-acid sequences. A falsy entry declines that row.
        v: Optional per-row V calls (same length). Each entry is an allele, comma-separated
            alleles, a list, or ``None`` to marginalize. ``None`` for the argument marginalizes
            every row.
        j: Optional per-row J calls (as ``v``).
        k: Distinct nucleotide candidates re-scored per row.
        threads: **Kernel threads** -- not worker processes. ``0`` = auto
            (``hardware_concurrency - 2``); batches under 64 rows stay single-threaded so their
            result is bitwise-identical to a serial run.
        resolve_genes: Resolve a **gene**-level call to a representative allele, as
            :func:`best_aa_scenarios`.

    Returns:
        Parallel columns, one entry per input row in input order: ``ok`` (0 for a row nothing
        explains), ``cdr3_nt``, ``v_call``, ``len_v``, ``j_call``, ``len_j``, ``d_call``,
        ``d_start``, ``d_end``, ``n_candidates``, ``pgen``, ``scenario_p``, ``runner_up_pgen``.
        The V/J/D names of a declined row are meaningless, not merely unset -- read ``ok`` first.

    Raises:
        KeyError: If any ``v``/``j`` entry names no gene the model carries.
        ValueError: If ``v`` or ``j`` is given with a length other than ``len(cdr3_aas)``.
    """
    import numpy as np

    from .._core import infer_nt_batch as _batch

    pm, vi, ji = pack(model)
    aas = [s.upper() if s else "" for s in cdr3_aas]
    if v is not None and len(v) != len(aas):
        raise ValueError("v must have the same length as cdr3_aas")
    if j is not None and len(j) != len(aas):
        raise ValueError("j must have the same length as cdr3_aas")
    alias = gene_to_allele(model) if resolve_genes else {}
    v_alts = [_alts(x, alias, vi, "V") for x in (v if v is not None else [None] * len(aas))]
    j_alts = [_alts(x, alias, ji, "J") for x in (j if j is not None else [None] * len(aas))]
    cols = _batch(pm, aas, v_alts, j_alts, k, threads)

    # Index -> allele name by fancy-indexing, so naming every row stays vectorized.
    names_v = np.array([a for a, _i in sorted(vi.items(), key=lambda kv: kv[1])], dtype=object)
    names_j = np.array([a for a, _i in sorted(ji.items(), key=lambda kv: kv[1])], dtype=object)
    names_d = np.array(model.genomic["genes_d"]["d_allele"].to_list()
                       if model.chain_type == "VDJ" else [], dtype=object)
    d = cols["d"]
    cols["v_call"] = names_v[np.maximum(cols.pop("v"), 0)]
    cols["j_call"] = names_j[np.maximum(cols.pop("j"), 0)]
    cols["d_call"] = (np.where(d >= 0, names_d[np.maximum(d, 0)], None) if names_d.size
                      else np.full(d.shape, None, dtype=object))
    del cols["d"]
    return cols
