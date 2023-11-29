import math


def test_get_volume(BaS_hybrid_calculation):
    assert math.isclose(BaS_hybrid_calculation.volume, 63.2552)


def test_get_energy(BaS_hybrid_calculation):
    assert math.isclose(BaS_hybrid_calculation.energy, -236035.1024)


def test_get_xc(BaS_hybrid_calculation):
    assert math.isclose(BaS_hybrid_calculation.xc, "hse06")

def test_get_num_atoms(BaS_hybrid_calculation):
    assert math.isclose(BaS_hybrid_calculation.num_atoms, 2)
