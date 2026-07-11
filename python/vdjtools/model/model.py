"""The :class:`Model` container — a manifest plus its polars marginal and germline tables."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl

from .schema import Manifest, validate_tables


@dataclass(slots=True)
class Model:
    """A V(D)J recombination model: declared graph + tidy polars tables.

    Args:
        manifest: Locus metadata and the recombination Bayes net.
        tables: Event name -> its long-format marginal ``pl.DataFrame``.
        genomic: ``"genes_v"`` / ``"genes_j"`` / (VDJ) ``"genes_d"`` -> germline reference frame.
    """

    manifest: Manifest
    tables: dict[str, pl.DataFrame]
    genomic: dict[str, pl.DataFrame]

    @property
    def locus(self) -> str:
        return self.manifest.locus

    @property
    def organism(self) -> str:
        return self.manifest.organism

    @property
    def chain_type(self) -> str:
        return self.manifest.chain_type

    def validate(self, *, tol: float = 1e-5) -> "Model":
        """Assert every event table has the right columns and normalizes; returns ``self``."""
        validate_tables(self.manifest, self.tables, tol=tol)
        return self

    def save(self, path: str | Path) -> None:
        """Write the model to a directory (``manifest.json`` + one parquet per table)."""
        from .io import save_model

        save_model(self, path)

    @classmethod
    def load(cls, path: str | Path) -> "Model":
        """Load a model previously written by :meth:`save`."""
        from .io import load_model

        return load_model(path)
