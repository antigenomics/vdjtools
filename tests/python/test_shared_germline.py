"""The germline API and bundled TRA models preserve shared alpha/delta alleles."""
import csv
from types import SimpleNamespace

import polars as pl
import pytest

from arda.paths import vdj_dir
from vdjtools.model import from_arda, load_bundled, load_germline


@pytest.mark.parametrize('organism', ['human', 'mouse'])
@pytest.mark.parametrize('legacy_anchors', [False, True])
def test_shared_germline_membership(organism, legacy_anchors, monkeypatch):
    if legacy_anchors:
        from arda import cdr3fix
        anchors = {key: SimpleNamespace(**{k: v for k, v in vars(anchor).items() if k != 'loci'})
                   for key, anchor in cdr3fix.load_anchors(organism).items()}
        monkeypatch.setattr(cdr3fix, 'load_anchors', lambda org: anchors)
    load_germline.cache_clear()
    with (vdj_dir(organism) / 'cdr3_anchors.tsv').open() as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    for locus in ('TRA', 'TRD'):
        expected = {(row['segment'], row['allele']) for row in rows if row['locus'] == locus}
        germline = load_germline(locus, organism)
        actual = set(germline.filter(pl.col('segment') != 'D').select('segment', 'allele').iter_rows())
        assert actual == expected, (locus, expected - actual, actual - expected)
    load_germline.cache_clear()


@pytest.mark.parametrize('organism,allele', [
    ('human', 'TRAV38-2/DV8*01'), ('mouse', 'TRAV13-4/DV7*01'),
])
def test_shared_allele_reaches_template_and_bundle(organism, allele):
    template = from_arda('TRA', organism)
    bundled = load_bundled('TRA', 'arda', organism=organism, collapse=False)
    for model in (template, bundled):
        assert allele in model.genomic['genes_v']['v_allele']
    assert set(bundled.genomic['genes_v']['v_allele']) == set(template.genomic['genes_v']['v_allele'])
    assert bundled.tables['v_choice'].filter(pl.col('v_allele') == allele)['p'].item() > 0
