from pathlib import Path
from runpy import run_path
import numpy as np
from quirksim.conformal_characteristics import ConformalString
from quirksim.decays import sample_decay
from quirksim.dynamics import C_CM_PER_PS as C


def test_flux_radial_matches_exact_constant_tension_solution():
    s = ConformalString.radial(225.,774.)
    y0=s.y.copy(); e0=np.hypot(225.,y0[11])
    assert s.advance(1.,max_step_ps=.05,integration_scheme='flux') == 'TARGET_REACHED'
    expected_p=y0[11]-s.force
    expected_x=y0[5]+C*(e0-np.hypot(225.,expected_p))/s.force
    np.testing.assert_allclose(s.y[11],expected_p,atol=1e-11,rtol=0)
    np.testing.assert_allclose(s.y[5],expected_x,atol=1e-13,rtol=0)
    assert s.snapshot()['chord_closure_error_cm'] < 1e-13


def test_flux_agrees_with_independent_exact_rotating_worldsheet():
    make=run_path(str(Path(__file__).resolve().parents[1]/'scripts/validate_conformal_boundary.py'))['rotating_string']
    s,omega,radius,_=make()
    assert s.advance(5.,max_step_ps=.1,integration_scheme='flux') == 'TARGET_REACHED'
    unit=np.array([np.cos(omega*5),np.sin(omega*5),0])
    np.testing.assert_allclose(s.y[2:8].reshape(2,3),[-radius*unit,radius*unit],atol=2e-12,rtol=0)


def test_flux_keeps_weak_jumps_causal_and_accounts_magnetic_impulse():
    s=ConformalString.boosted_radial(225.2,774.)
    for endpoint in (0,1):
        products=sample_decay('pion',s.y[8+3*endpoint:11+3*endpoint],225.2,225.,np.random.default_rng(endpoint))
        s.apply_endpoint_decay(endpoint,225.,products)
    initial=s.snapshot()
    assert s.advance(5.,max_step_ps=.025,causal_fraction=.01,integration_scheme='flux') == 'TARGET_REACHED'
    final=s.snapshot()
    assert s.reflections and s.emissions
    assert abs(final['total_energy_gev']-initial['total_energy_gev']) < 1e-9
    np.testing.assert_allclose(np.array(final['total_momentum_gev'])-final['external_impulse_gev'],
                              np.array(initial['total_momentum_gev'])-initial['external_impulse_gev'],atol=1e-9,rtol=0)
    assert final['chord_closure_error_cm'] < 2e-12


def test_flux_matching_surface_stops_without_annihilation():
    s=ConformalString.radial(225.,774.,length_cm=.001,beta=0.)
    assert s.advance(10.,max_step_ps=.1,matching_radius_cm=.0005,
                     integration_scheme='flux') == 'QUANTUM_MATCHING_SURFACE'
    assert abs(np.linalg.norm(s.y[5:8]-s.y[2:5])-.0005) < 2e-12
    assert s.emissions == []
