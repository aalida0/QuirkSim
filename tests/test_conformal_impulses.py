"""Local weak matching, magnetic work and finite signal propagation."""
import numpy as np
import pytest
from quirksim.conformal_characteristics import ConformalString, UnitRayHistory
from quirksim.physics import benchmark
from quirksim.decays import sample_decay
from quirksim.conformal_evolution import evolve


def test_one_sided_null_jump_has_zero_integrated_impulse():
    rays = UnitRayHistory([0., 1.], [[1., 0, 0]]*2)
    rays.jump(1., [0., 1., 0.])
    rays.append(2., [0., 1., 0.])
    np.testing.assert_allclose(rays.value(1-1e-8), [1., 0, 0])
    np.testing.assert_allclose(rays.value(1.), [0., 1., 0])
    np.testing.assert_allclose(rays.primitive(1.), [1., 0, 0])
    np.testing.assert_allclose(rays.primitive(2.), [1., 1., 0])


@pytest.mark.parametrize('channel', ['pion', 'electron_neutrino', 'muon_neutrino'])
@pytest.mark.parametrize('endpoint', [0, 1])
def test_local_weak_decay_conserves_four_momentum_and_leaves_bulk_integrals(channel, endpoint):
    m = benchmark()['charged_mass_gev']
    s = ConformalString.boosted_radial(m, 774., b_tesla=2.)
    s.advance(.5, .05)
    before = s.snapshot()
    p = s.y[8+3*endpoint:11+3*endpoint].copy()
    products = sample_decay(channel, p, m, 225., np.random.default_rng(96))
    record = s.apply_endpoint_decay(endpoint, 225., products)
    after = s.snapshot()
    emitted = sum((v for k, v in products.items() if k != 'daughter'), np.zeros(4))
    assert s.charges[endpoint] == 0.
    np.testing.assert_allclose(np.r_[before['total_energy_gev'], before['total_momentum_gev']],
                              np.r_[after['total_energy_gev'], after['total_momentum_gev']]+emitted, atol=1e-10)
    assert before['string_energy_gev'] == after['string_energy_gev']
    assert max(np.abs(record['four_momentum_residual_gev'])) < 1e-10


def test_recoil_signal_cannot_change_the_other_endpoint_before_light_travel():
    m = benchmark()['charged_mass_gev']
    kicked = ConformalString.radial(m, 774., length_cm=.3, beta=0.)
    control = ConformalString.radial(m, 774., length_cm=.3, beta=0.)
    products = sample_decay('pion', [0., 0., 0.], m, 225., np.random.default_rng(12))
    kicked.apply_endpoint_decay(0, 225., products)
    kicked.advance(2., .02)
    control.advance(2., .02)
    np.testing.assert_allclose(kicked.y[5:8], control.y[5:8], atol=1e-15)
    np.testing.assert_allclose(kicked.y[11:14], control.y[11:14], atol=1e-14)
    assert abs(kicked.y[9])+abs(kicked.y[10]) > 1e-5
    assert kicked.advance(30., .02) == 'TARGET_REACHED'
    assert len(kicked.reflections) >= 3
    assert kicked.snapshot()['chord_closure_error_cm'] < 1e-10


def test_magnetic_work_is_zero_and_momentum_impulse_is_accounted():
    s = ConformalString.boosted_radial(225., 774.)
    initial = s.snapshot()
    s.advance(5., .05)
    z = s.snapshot()
    assert abs(z['total_energy_gev']-initial['total_energy_gev']) < 1e-9
    np.testing.assert_allclose(np.array(z['total_momentum_gev'])-z['external_impulse_gev'],
                              initial['total_momentum_gev'], atol=1e-8)


def test_boosted_early_weak_pilot_reflects_one_sided_wave_without_smeared_impulse():
    # This catches using the right-hand incident ray in the final RK stage at
    # an incoming jump, which had generated a false ~3e-8 cm geometry jump.
    r = evolve(seed=25, duration_ps=30., max_step_ps=.05,
               max_ray_angle=.0015625, causal_fraction=.003125)
    assert r['termination'] == 'TARGET_REACHED'
    assert len(r['decays']) == 2
    assert len(r['wave_reflections']) >= 2
    assert max(z['chord_closure_error_cm'] for z in r['records']) < 1e-10
    assert max(abs(x) for x in r['momentum_with_weak_products_minus_field_residual_gev']) < 1e-9
    assert all(d['parent_charge'] == d['daughter_charge']+d['emitted_charge'] for d in r['decays'])
