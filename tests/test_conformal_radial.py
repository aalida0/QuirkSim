import numpy as np
import pytest
from quirksim.conformal_radial import signed_radial_boundary
from quirksim.radial_ng_reference import ExactRadialNG


def test_four_crossings_keep_endpoint_identity_momentum_and_proper_time():
    ref = ExactRadialNG()
    segments, crossings = signed_radial_boundary(ref.mass, 774., .1, .4, 1000., 2.)
    assert len(crossings) == 4
    for i, crossing in enumerate(crossings):
        assert crossing['time_ps'] == pytest.approx(
            ref.first_crossing_ps+i*ref.crossing_interval_ps, abs=1e-9)
        np.testing.assert_allclose(segments[i].y[:, -1], segments[i+1].y[:, 0], atol=1e-12)
    r, p, tau = segments[-1].y[:, -1]
    exact = ref.state(1000.)
    np.testing.assert_allclose([r, p, tau], [exact['right_position_cm'],
        exact['right_momentum_gev'], exact['proper_time_each_ps']], atol=1e-9)
