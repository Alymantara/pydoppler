import numpy as np

from pydoppler import pydoppler as pydoppler_module


def test_stream_calculate_preserves_mass_ratio_near_unity():
    exact = pydoppler_module.stream_calculate(1.0, ni=30, nj=30)
    near = pydoppler_module.stream_calculate(0.9999, ni=30, nj=30)
    assert not np.array_equal(near[0], exact[0])
    assert np.allclose(near[0], exact[0], rtol=0, atol=3e-5)
    assert pydoppler_module.rlq1(1.0) == 0.5
