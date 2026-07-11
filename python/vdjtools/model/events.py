"""V(D)J recombination events — the Bayes-net structure, declared as data.

A generative model's factorization is a set of named *events*, each conditioned on zero or
more parent events (its ``given`` list). This mirrors IGoR's ``@Edges`` and OLGA's parameter
factorization, but is stored declaratively in the model manifest rather than baked into a
code path — so a locus can change its conditioning (say ``P(D | J)`` vs ``P(D | V, J)``, or
add a tandem-D slot) without touching the engine.

This module is just the vocabulary and validation. The concrete graph for a given model
lives in its :class:`~vdjtools.model.schema.Manifest` (populated by the loader or by EM).
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class EventKind(str, Enum):
    """The kind of a recombination event — fixes the realization columns of its table.

    - ``GENE_CHOICE``   — pick a germline allele (V, D, or J).            realization: ``<seg>_allele``
    - ``N_D``           — number of D segments in the junction (0/1/2).   realization: ``n_d``
    - ``DELETION``      — 3'/5' exonuclease trimming of one segment end.  realization: ``ndel`` (may be <0: P-nt)
    - ``DELETION_2D``   — joint 5'+3' trimming of a D segment.            realization: ``ndel5, ndel3``
    - ``INS_LENGTH``    — number of non-templated (N-region) insertions.  realization: ``length``
    - ``DINUCLEOTIDE``  — Markov transition of the N-region nt sequence.  realization: ``from_nt, to_nt``
    """

    GENE_CHOICE = "gene_choice"
    N_D = "n_d"
    DELETION = "deletion"
    DELETION_2D = "deletion_2d"
    INS_LENGTH = "ins_length"
    DINUCLEOTIDE = "dinucleotide"


@dataclass(frozen=True, slots=True)
class Event:
    """One node of the recombination Bayes net.

    Args:
        name: Unique event name (also the marginal table's stem, e.g. ``"v_choice"``).
        kind: The :class:`EventKind`, which fixes the table's realization columns.
        given: Names of parent events this one is conditioned on (its normalization key).
    """

    name: str
    kind: EventKind
    given: tuple[str, ...] = ()


def validate_graph(events: dict[str, Event]) -> None:
    """Check that an event graph is well-formed: parents exist and there are no cycles.

    Args:
        events: Mapping of event name to :class:`Event`.

    Raises:
        ValueError: If a ``given`` names a missing event, or the graph has a cycle.
    """
    for ev in events.values():
        for parent in ev.given:
            if parent not in events:
                raise ValueError(f"event {ev.name!r} is given unknown parent {parent!r}")

    # DFS cycle check over the `given` edges.
    WHITE, GREY, BLACK = 0, 1, 2
    color = {name: WHITE for name in events}

    def visit(name: str, stack: tuple[str, ...]) -> None:
        if color[name] == GREY:
            cycle = " -> ".join((*stack, name))
            raise ValueError(f"event graph has a cycle: {cycle}")
        if color[name] == BLACK:
            return
        color[name] = GREY
        for parent in events[name].given:
            visit(parent, (*stack, name))
        color[name] = BLACK

    for name in events:
        visit(name, ())
