import numpy as np
import pytest
from vdjtools.model.selection import ALICE_TRB_SELECTION_Q, selection_mass

def test_scaling_and_validation():
    np.testing.assert_array_equal(selection_mass([0.,.1],4),[0.,.4])
    assert selection_mass(.1)==.1
    assert ALICE_TRB_SELECTION_Q==9.41
    for mass,q in [(-1,1),(float('nan'),1),(.1,0),(.1,float('inf'))]:
        with pytest.raises(ValueError):selection_mass(mass,q)
