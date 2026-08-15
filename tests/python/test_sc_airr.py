"""Tests for the single-cell AIRR Rearrangement interchange layer (vdjtools.sc.airr).

This is the contract every downstream bridge (scirpy, dandelion, scRepertoire) is built
on, so it is pinned with NO optional dependencies -- these tests must never skip.

Pins: the emitted column set and AIRR spelling (junction, not junction_nt), sequence_id
synthesis for frames that lack one, that consensus_count carries the read count that
scRepertoire's parser reads, the to_airr/from_airr round-trip, and the scRepertoire
export shapes. Plus the io.sniff_format regression: a barcoded AIRR table must not be
read as bulk, because that silently discards cell_id.
"""
import polars as pl
import pytest

from vdjtools import sc
from vdjtools.sc.airr import AIRR_COLUMNS, CONSENSUS_COUNT, JUNCTION
from vdjtools.sc.read import SC_COLUMNS


def _cells():
    """A minimal sc long frame: two cells, one with two chains."""
    return pl.DataFrame({
        "cell_id": ["c1", "c1", "c2"],
        "sequence_id": ["c1_contig_1", "c1_contig_2", "c2_contig_1"],
        "locus": ["TRA", "TRB", "TRB"],
        "v_call": ["TRAV1-2", "TRBV20-1", "TRBV9"],
        "d_call": [None, "TRBD1", None],
        "j_call": ["TRAJ33", "TRBJ2-7", "TRBJ2-3"],
        "c_call": ["TRAC", "TRBC2", "TRBC2"],
        "junction_aa": ["CAVRDSNYQLIW", "CASSLGQAYEQYF", "CASSVDTQYF"],
        "junction_nt": ["TGTGCC", "TGTGCCAGC", "TGTGCCAGCAGC"],
        "duplicate_count": [120, 340, 88],
        "umi_count": [6, 14, 4],
        "clone_id": ["clonotype1", "clonotype1", "clonotype2"],
        "productive": [True, True, True],
    }).select(SC_COLUMNS)


def test_to_airr_emits_the_airr_column_set_and_spelling():
    airr = sc.to_airr(_cells())
    assert airr.columns == AIRR_COLUMNS
    # AIRR spells the nucleotide junction `junction`; vdjtools stores it as junction_nt.
    assert JUNCTION in airr.columns and "junction_nt" not in airr.columns
    assert airr[JUNCTION].to_list() == ["TGTGCC", "TGTGCCAGC", "TGTGCCAGCAGC"]
    assert airr.height == 3


def test_consensus_count_carries_the_read_count():
    # scRepertoire's .parseAIRR reads consensus_count as `reads`; scirpy/dandelion prefer
    # umi_count. Emitting both is what lets one file feed all three.
    airr = sc.to_airr(_cells())
    assert airr[CONSENSUS_COUNT].to_list() == [120, 340, 88]
    assert airr["duplicate_count"].to_list() == [120, 340, 88]
    assert airr["umi_count"].to_list() == [6, 14, 4]


def test_sequence_id_is_synthesised_per_cell_when_absent():
    # dandelion derives cell_id back out of `<cell>_contig_<n>`, so the round-trip closes
    # even for a frame that never carried a sequence_id.
    cells = _cells().drop("sequence_id")
    airr = sc.to_airr(cells)
    assert airr["sequence_id"].to_list() == ["c1_contig_1", "c1_contig_2", "c2_contig_1"]


def test_sequence_id_synthesis_fills_only_the_nulls():
    cells = _cells().with_columns(
        pl.Series("sequence_id", ["given", None, None])
    )
    airr = sc.to_airr(cells)
    assert airr["sequence_id"].to_list() == ["given", "c1_contig_2", "c2_contig_1"]


def test_missing_optional_columns_become_nulls_not_errors():
    cells = _cells().drop("d_call", "clone_id", "umi_count")
    airr = sc.to_airr(cells)
    assert airr.columns == AIRR_COLUMNS
    assert airr["d_call"].null_count() == 3
    assert airr["umi_count"].null_count() == 3


def test_to_airr_requires_cell_id():
    with pytest.raises(ValueError, match="cell_id"):
        sc.to_airr(_cells().drop("cell_id"))


def test_from_airr_round_trips_the_canonical_frame():
    cells = _cells()
    assert sc.from_airr(sc.to_airr(cells)).equals(cells)


def test_from_airr_accepts_either_junction_spelling():
    cells = _cells()
    airr = sc.to_airr(cells)
    # A producer that used the vdjtools spelling is still readable.
    renamed = airr.rename({JUNCTION: "junction_nt"})
    assert sc.from_airr(renamed)["junction_nt"].to_list() == cells["junction_nt"].to_list()


def test_from_airr_falls_back_to_consensus_count_for_reads():
    airr = sc.to_airr(_cells()).drop("duplicate_count")
    assert sc.from_airr(airr)["duplicate_count"].to_list() == [120, 340, 88]


def test_write_airr_is_readable_back_through_read_airr_cell(tmp_path):
    out = sc.write_airr(_cells(), tmp_path / "airr_rearrangement.tsv")
    back = sc.read_airr_cell(out)
    assert back.columns == SC_COLUMNS
    assert back["cell_id"].to_list() == ["c1", "c1", "c2"]
    assert back["junction_nt"].to_list() == ["TGTGCC", "TGTGCCAGC", "TGTGCCAGCAGC"]
    assert back["duplicate_count"].to_list() == [120, 340, 88]


def test_write_screpertoire_airr_has_exactly_the_columns_parseairr_wants(tmp_path):
    out = sc.write_screpertoire(_cells(), tmp_path / "airr_rearrangement.tsv")
    header = pl.read_csv(out, separator="\t", n_rows=1).columns
    assert header == ["cell_id", "locus", "consensus_count", "v_call", "d_call",
                      "j_call", "c_call", "junction", "junction_aa"]


def test_write_screpertoire_10x_has_the_contig_annotation_columns(tmp_path):
    out = sc.write_screpertoire(_cells(), tmp_path / "contigs.csv", format="10x")
    header = pl.read_csv(out, n_rows=1).columns
    assert header == ["barcode", "contig_id", "chain", "v_gene", "d_gene", "j_gene",
                      "c_gene", "cdr3", "cdr3_nt", "reads", "umis", "productive"]


def test_write_screpertoire_rejects_an_unknown_format(tmp_path):
    with pytest.raises(ValueError, match="airr.*10x"):
        sc.write_screpertoire(_cells(), tmp_path / "x.tsv", format="seurat")


def test_barcoded_airr_is_not_silently_read_as_bulk(tmp_path):
    """Regression: io.read used to collapse a single-cell table, dropping cell_id."""
    from vdjtools import io

    path = sc.write_airr(_cells(), tmp_path / "airr_rearrangement.tsv")
    assert io.sniff_format(path) == "airr_cell"
    with pytest.raises(ValueError, match="single-cell"):
        io.read(path)
    # Pooling it into a bulk repertoire on purpose is still allowed.
    pooled = io.read(path, fmt="airr")
    assert "cell_id" not in pooled.columns and pooled.height > 0
