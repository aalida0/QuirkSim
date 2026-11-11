import numpy as np
import pytest
from quirksim.conformal_characteristics import ConformalString, UnitRayHistory
from quirksim.radial_ng_reference import ExactRadialNG


def pulse_string(samples=121):
    s = np.linspace(0., 6., samples)
    angle = np.where((s > 2.) & (s < 3.), .03*np.sin(np.pi*(s-2))**4, 0.)
    a = np.c_[np.cos(angle), np.sin(angle), np.zeros_like(s)]
    b = np.tile([-1., 0, 0], (len(s), 1))
    return ConformalString([225., 225.], 774., s, a, b, np.zeros((2, 3)))


def test_history_null_interpolation_and_exact_integral():
    h = UnitRayHistory([0., 1.], [[1., 0, 0], [0., 1., 0]])
    np.testing.assert_allclose(h.value(.5), [2**-.5, 2**-.5, 0], atol=1e-15)
    np.testing.assert_allclose(h.primitive(1), [2/np.pi, 2/np.pi, 0], atol=1e-15)
    with pytest.raises(ValueError, match='causal'):
        h.value(1.01)


def test_characteristic_radial_matches_independent_exact_reference_before_contact():
    ref = ExactRadialNG()
    solver = ConformalString.radial(ref.mass, 774.)
    for t in (10., 50., 100., 150.):
        assert solver.advance(t, .5) == 'TARGET_REACHED'
        z = solver.snapshot()
        exact = ref.state(t)
        assert z['positions_cm'][1][0] == pytest.approx(exact['right_position_cm'], abs=2e-11)
        assert z['momenta_gev'][1][0] == pytest.approx(exact['right_momentum_gev'], abs=1e-10)
        assert z['proper_times_ps'][1] == pytest.approx(exact['proper_time_each_ps'], abs=2e-9)
        assert z['total_energy_gev'] == pytest.approx(ref.energy_total, abs=2e-9)
        assert z['chord_closure_error_cm'] < 1e-12


def test_pulse_hits_massive_end_then_reflects_causally_to_other_end():
    solver = pulse_string()
    initial = solver.snapshot()
    solver.advance(1.8, .01)
    assert abs(solver.y[9]) < 1e-15
    solver.advance(3.5, .01)
    assert solver.y[9] > .005
    solver.advance(6., .01)
    assert abs(solver.y[12]) < 1e-15
    solver.advance(9., .01)
    assert abs(solver.y[12]) > .001
    z = solver.snapshot()
    assert z['max_conformal_norm_residual'] < 2e-15
    assert z['max_conformal_orthogonality_residual'] < 2e-15
    assert abs(z['total_energy_gev']-initial['total_energy_gev']) < 2e-6
    np.testing.assert_allclose(z['total_momentum_gev'], initial['total_momentum_gev'], atol=5e-6)
    assert z['chord_closure_error_cm'] < 5e-7


def test_chart_collapse_stops_without_endpoint_relabeling():
    solver = ConformalString.radial(225., 774., length_cm=1e-5, beta=0.)
    status = solver.advance(10., .02, minimum_interval_ps=1e-7)
    assert status == 'STOP_ENERGY_COORDINATE_COLLAPSE'
    assert solver.y[1] > solver.y[0]
    assert solver.y[5] > solver.y[2]
    assert solver.t < 10.


def test_incompatible_initial_velocity_and_invalid_step_are_rejected():
    with pytest.raises(ValueError, match='incompatible'):
        ConformalString([225., 225.], 774., [0, 1], [[1., 0, 0]]*2,
                        [[-1., 0, 0]]*2, [[0., 2., 0], [0., 0, 0]])
    with pytest.raises(ValueError):
        pulse_string().advance(1., 0.)
