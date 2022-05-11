import math

def test_get_volume(BaS_hybrid_calculation):
	assert math.isclose(BaS_hybrid_calculation.volume,63.2552)
