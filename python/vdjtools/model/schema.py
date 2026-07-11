"""Polars schema for a V(D)J recombination model + its manifest.

A model is a directory: one tidy (long-format) parquet per event marginal, parquet tables
for the germline references, and a ``manifest.json`` that declares the Bayes-net graph
(:mod:`~vdjtools.model.events`) plus locus metadata. This is the clean tabular replacement
for IGoR's ``model_parms.txt`` / ``model_marginals.txt`` grammar: every probability row is
self-describing and every M-step normalization is one ``group_by(key).over(...)``.

Conventions (match OLGA so a loaded model is bit-faithful — see the loader):

- Nucleotides are integer-coded ``A,C,G,T = 0,1,2,3`` everywhere.
- ``ndel`` is the **biological** deletion count: negative values are palindromic (P-) nt,
  ``0`` is a flush cut, positive values trim germline. (OLGA stores ``ndel + max_palindrome``
  as an array index; we store the biological value.)
- A dinucleotide table row ``(from_nt, to_nt, p)`` is ``p = P(next = to_nt | prev = from_nt)``,
  so it normalizes within ``from_nt`` (OLGA's column-stochastic ``R[next, prev]``).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

import polars as pl

from .events import Event, EventKind, validate_graph

# Realization columns (with dtypes) contributed by each event kind, in table order.
# The event's ``given`` columns and the ``p`` column are appended by :func:`table_columns`.
_REALIZATION: dict[EventKind, dict[str, pl.DataType]] = {
    EventKind.GENE_CHOICE: {},  # the chosen-allele column is named per event (see _allele_col)
    EventKind.N_D: {"n_d": pl.UInt8},
    EventKind.DELETION: {"ndel": pl.Int16},
    EventKind.DELETION_2D: {"ndel5": pl.Int16, "ndel3": pl.Int16},
    EventKind.INS_LENGTH: {"length": pl.Int16},
    EventKind.DINUCLEOTIDE: {"from_nt": pl.UInt8, "to_nt": pl.UInt8},
}


def _allele_col(event: Event) -> str:
    """The chosen-allele column name for a gene-choice event (``v_choice`` -> ``v_allele``)."""
    seg = event.name.split("_")[0]  # v_choice->v, j_choice->j, d_gene->d, d2_gene->d2
    return f"{seg}_allele"


def table_columns(event: Event) -> dict[str, pl.DataType]:
    """Full column schema (name -> dtype) for an event's marginal table.

    Layout: ``given`` allele columns, then the event's own realization columns, then ``p``.
    """
    cols: dict[str, pl.DataType] = {g: pl.Utf8 for g in _given_allele_cols(event)}
    if event.kind is EventKind.GENE_CHOICE:
        cols[_allele_col(event)] = pl.Utf8
    else:
        cols.update(_REALIZATION[event.kind])
    cols["p"] = pl.Float64
    return cols


def _given_allele_cols(event: Event) -> list[str]:
    """Column names for the event's parents (each parent is a gene-choice → its allele col)."""
    return [f"{g.split('_')[0]}_allele" for g in event.given]


def normalization_keys(event: Event) -> list[str]:
    """Columns the table's ``p`` must sum to 1 within (the M-step / validation group key).

    Parents (``given``) always; plus ``from_nt`` for dinucleotide tables (column-stochastic).
    """
    keys = _given_allele_cols(event)
    if event.kind is EventKind.DINUCLEOTIDE:
        keys = [*keys, "from_nt"]
    return keys


@dataclass(frozen=True, slots=True)
class Manifest:
    """Model metadata + the declared recombination Bayes net.

    Args:
        locus: e.g. ``"TRB"``.
        organism: e.g. ``"human"``.
        chain_type: ``"VDJ"`` (has D) or ``"VJ"`` (no D).
        events: The recombination graph, name -> :class:`~vdjtools.model.events.Event`.
        palindrome_max: Max palindromic nt per trimmable end (e.g. ``{"v_3": 4, "j_5": 4}``).
        model_version: Schema/model version tag.
        source: Free-text provenance (e.g. ``"olga:human_T_beta"``).
        error_rate: Optional per-nt error rate (unused by Pgen; carried for round-trip).
    """

    locus: str
    organism: str
    chain_type: str
    events: dict[str, Event]
    palindrome_max: dict[str, int] = field(default_factory=dict)
    model_version: str = "2.0.0"
    source: str = ""
    error_rate: float | None = None

    def __post_init__(self) -> None:
        if self.chain_type not in ("VDJ", "VJ"):
            raise ValueError(f"chain_type must be 'VDJ' or 'VJ', got {self.chain_type!r}")
        validate_graph(self.events)

    def to_json(self) -> str:
        obj = {
            "locus": self.locus,
            "organism": self.organism,
            "chain_type": self.chain_type,
            "model_version": self.model_version,
            "source": self.source,
            "error_rate": self.error_rate,
            "palindrome_max": self.palindrome_max,
            "events": {
                name: {"kind": ev.kind.value, "given": list(ev.given)}
                for name, ev in self.events.items()
            },
        }
        return json.dumps(obj, indent=2)

    @classmethod
    def from_json(cls, text: str) -> "Manifest":
        obj = json.loads(text)
        events = {
            name: Event(name=name, kind=EventKind(spec["kind"]), given=tuple(spec["given"]))
            for name, spec in obj["events"].items()
        }
        return cls(
            locus=obj["locus"],
            organism=obj["organism"],
            chain_type=obj["chain_type"],
            events=events,
            palindrome_max=obj.get("palindrome_max", {}),
            model_version=obj.get("model_version", "2.0.0"),
            source=obj.get("source", ""),
            error_rate=obj.get("error_rate"),
        )


def validate_tables(manifest: Manifest, tables: dict[str, pl.DataFrame], *, tol: float = 1e-5) -> None:
    """Check every event table has the right columns and normalizes within its key.

    Args:
        manifest: The model manifest declaring the events.
        tables: Event name -> its marginal ``pl.DataFrame``.
        tol: Absolute tolerance for the "sums to 1" check.

    Raises:
        ValueError: On a missing table, wrong columns, or a group whose ``p`` ≠ 1.
    """
    for name, event in manifest.events.items():
        if name not in tables:
            raise ValueError(f"model is missing the marginal table for event {name!r}")
        df = tables[name]
        expected = table_columns(event)
        if set(df.columns) != set(expected):
            raise ValueError(
                f"table {name!r} columns {sorted(df.columns)} != expected {sorted(expected)}"
            )
        keys = normalization_keys(event)
        sums = (
            df.group_by(keys).agg(pl.col("p").sum().alias("s"))
            if keys
            else df.select(pl.col("p").sum().alias("s"))
        )
        # Each group's p sums to 1, or to 0 for an *undefined conditional* — a parent value
        # (gene/nt) that never occurs, so its conditional is empty. OLGA keeps such all-zero
        # columns for gene-index alignment; a faithful load preserves them.
        bad = sums.filter((pl.col("s") - 1.0).abs() > tol).filter(pl.col("s").abs() > tol)
        if bad.height:
            raise ValueError(
                f"table {name!r} has {bad.height} group(s) whose probabilities sum to neither "
                f"1 nor 0 (within {tol}); first offender: {bad.row(0)}"
            )
