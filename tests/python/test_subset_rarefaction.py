"""Subset discovery uses total reads, with exact finite-population oracles."""
import math
import numpy as np
import pytest
from vdjtools.stats import subset_rarefaction_batch, gini_coefficient


def test_subset_batch_matches_exact_sampling_and_endpoints():
    vectors=[[1,2],[],[8],[1]*8,[2,3]]
    curve=subset_rarefaction_batch(vectors,[8]*5,sizes=range(9),threads=2)
    assert curve.height==45
    for i,x in enumerate(vectors):
        got=curve.filter(curve['sample']==i)['qD'].to_numpy()
        want=[sum(1-(math.comb(8-c,m) if 8-c>=m else 0)/math.comb(8,m) for c in x) for m in range(9)]
        np.testing.assert_allclose(got,want,atol=1e-12)
    assert gini_coefficient(np.array([1,1,1]))==0
    assert gini_coefficient(np.array([1,1,4]))==pytest.approx(1/3)
    assert gini_coefficient(np.array([]))==0
    assert gini_coefficient(np.array([4]))==0


@pytest.mark.parametrize('counts,total,depths', [([1.5],8,[1]),([9],8,[1]),([1],0,[0]),([1],8,[9]),([1],8,[.5]),([float('nan')],8,[1])])
def test_subset_rejects_invalid_inputs(counts,total,depths):
    with pytest.raises(ValueError):
        subset_rarefaction_batch([counts],[total],sizes=depths,threads=1)
