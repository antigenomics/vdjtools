"""Pool clonotypes across samples (pure polars).

Reimplements the legacy ``PoolSamples`` / ``SampleAggregator`` / ``PooledSample``:
collapse clonotypes across a set of samples by a chosen match key, summing counts
and recomputing frequency, and annotate each pooled clonotype with its incidence,
occurrence count and convergence.
"""
from __future__ import annotations

import polars as pl

from ..io.schema import (
    JUNCTION_AA,
    JUNCTION_NT,
    COUNT,
    D_CALL,
    C_CALL,
    J_CALL,
    V_CALL,
    normalize,
)

#: Legacy ``OverlapType`` match keys -> the clonotype columns they compare on.
KEY_SETS: dict[str, tuple[str, ...]] = {
    "strict": (JUNCTION_NT, V_CALL, J_CALL),   # OverlapType.Strict
    "nt": (JUNCTION_NT,),                       # Nucleotide
    "ntV": (JUNCTION_NT, V_CALL),               # NucleotideV
    "ntVJ": (JUNCTION_NT, V_CALL, J_CALL),      # NucleotideVJ
    "aa": (JUNCTION_AA,),                       # AminoAcid
    "aaV": (JUNCTION_AA, V_CALL),               # AminoAcidV
    "aaVJ": (JUNCTION_AA, V_CALL, J_CALL),      # AminoAcidVJ
}


def resolve_key(key: "str | tuple[str, ...]") -> list[str]:
    """Resolve a match-key name (or explicit column tuple) to key columns.

    Args:
        key: One of the legacy overlap-type names (``"strict"``, ``"nt"``,
            ``"ntV"``, ``"ntVJ"``, ``"aa"``, ``"aaV"``, ``"aaVJ"``) or an explicit
            tuple/list of column names.

    Returns:
        The list of clonotype columns forming the key.

    Raises:
        ValueError: If ``key`` is an unknown name.
    """
    if isinstance(key, str):
        if key not in KEY_SETS:
            raise ValueError(f"unknown key {key!r}; expected one of {sorted(KEY_SETS)}")
        return list(KEY_SETS[key])
    return list(key)


def pool_samples(samples: "list[pl.DataFrame]", key: "str | tuple[str, ...]" = "aa",
                 sample_col: str | None = None) -> pl.DataFrame:
    """Pool clonotypes across samples, summing counts and recomputing frequency.

    Reimplements ``PooledSample`` (legacy default match key ``"aa"``). For each
    distinct match key across all samples the pooled clonotype carries:

    - ``duplicate_count`` — summed read count across all samples.
    - ``frequency`` — pooled count over the pool's total reads.
    - ``incidence`` — number of distinct samples the clonotype occurs in.
    - ``occurrences`` — total number of clonotype rows aggregated (a clonotype can
      appear once per sample, and, for amino-acid keys, via several nucleotide
      variants within one sample).
    - ``convergence`` — number of distinct nucleotide variants (``junction_nt`` +
      ``v_call`` + ``j_call``, the legacy strict key) collapsed into the pooled
      clonotype. This is 1 for nucleotide-level keys and counts convergent
      recombination for amino-acid-level keys.

    The representative non-key fields (e.g. the ``junction_nt`` of an amino-acid pool)
    are taken from the most abundant contributing row (legacy
    ``MaxClonotypeAggregator``).

    Args:
        samples: A list of per-sample clonotype frames, **or** a single long frame
            (pass ``[df]``) split by ``sample_col``.
        key: Match key name or explicit column tuple (see :func:`resolve_key`).
        sample_col: If the input is a single long frame carrying a sample-id
            column, its name; used only to count ``incidence`` correctly.

    Returns:
        A pooled clonotype frame in the canonical schema plus ``incidence``,
        ``occurrences`` and ``convergence``, sorted by descending
        ``duplicate_count``.
    """
    key = resolve_key(key)
    frames = []
    for i, s in enumerate(samples):
        sid = (s[sample_col].to_list() if (sample_col and sample_col in s.columns)
               else [i] * s.height)
        frames.append(normalize(s).with_columns(pl.Series("_sample", sid)))
    long = pl.concat(frames, how="vertical_relaxed")

    fields = (V_CALL, D_CALL, J_CALL, C_CALL, JUNCTION_AA, JUNCTION_NT)
    strict = [JUNCTION_NT, V_CALL, J_CALL]
    agg = (
        long.group_by(key, maintain_order=True)
        .agg(
            pl.col(COUNT).sum().alias(COUNT),
            pl.col("_sample").n_unique().alias("incidence"),
            pl.len().alias("occurrences"),
            pl.struct(strict).n_unique().alias("convergence"),
            # representative: fields of the single most-abundant contributing row
            *[pl.col(c).sort_by(COUNT, descending=True).first().alias(f"_rep_{c}")
              for c in fields],
        )
    )
    # Fill each canonical field from the representative where it is not the key.
    agg = agg.with_columns(
        [(pl.col(c) if c in key else pl.col(f"_rep_{c}")).alias(c) for c in fields]
    )

    out = normalize(agg, recompute_freq=True)
    out = out.with_columns(agg["incidence"], agg["occurrences"], agg["convergence"])
    return out.sort(COUNT, descending=True, maintain_order=True)
