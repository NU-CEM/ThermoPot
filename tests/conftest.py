import pytest
import os
from thermopot installcalculations, materials

@pytest.fixture()
def BaS_hybrid_calculation():
    return calculations.calculation(os.path.join(os.path.dirname(
        __file__), 'data/BaS_hse06_outfile'))

@pytest.fixture()
def BaS_pbesol_calculation():
    return calculations.calculation(os.path.join(os.path.dirname(
        __file__), 'data/BaS_pbesol_relax_outfile'))

@pytest.fixture()
def BaS_solid():
    return materials.solid('BaS',{'Ba':1, 'S':1},os.path.join(os.path.dirname(
        __file__), 'data/BaS_phonons'), BaS_pbesol_calculation)

@pytest.fixture()
def BaS_BaS_reaction():
    return reactions.reaction({BaS_solid: 1, BaS_solid : 1})

