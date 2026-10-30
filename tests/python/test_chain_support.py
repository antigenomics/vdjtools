import polars as pl
import pytest
from vdjtools.sc import chain_support

def test_pair_target_union_intersection_and_duplicate_hits():
    hits=pl.DataFrame(dict(receptor=['p','p','p','p','orphan'],target=['X','X','X','Y','X'],locus=['TRA','TRA','TRB','TRA','TRB']))
    result=chain_support(hits,by=['receptor','target']).sort('receptor','target')
    assert result.select('tra_supported','trb_supported','union_supported','intersection_supported').rows()==[(False,True,True,False),(True,True,True,True),(True,False,True,False)]
    assert chain_support(hits.head(0),by=['receptor','target']).height==0
    with pytest.raises(ValueError):chain_support(hits.with_columns(pl.lit('IGH').alias('locus')),by=['receptor','target'])
