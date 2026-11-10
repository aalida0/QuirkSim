import numpy as np
import pytest
from pathlib import Path
from runpy import run_path
from quirksim.conformal_characteristics import ConformalString
from quirksim.dynamics import C_CM_PER_PS as C

rotating_string = run_path(str(Path(__file__).resolve().parents[1]/
                              'scripts/validate_conformal_boundary.py'))['rotating_string']


def test_quantum_matching_censors_before_the_classical_chart_collapses():
    solver = ConformalString.radial(225., 774., length_cm=.001, beta=0.)
    status = solver.advance(10., .1, matching_radius_cm=.0005)
    assert status == 'QUANTUM_MATCHING_SURFACE'
    assert np.linalg.norm(solver.y[5:8]-solver.y[2:5]) == pytest.approx(.0005, abs=1e-12)
    assert solver.y[1] > solver.y[0]
    assert solver.emissions == []


def test_classical_poincare_charges_match_independent_rotating_solution():
    solver, omega, radius, momentum = rotating_string()
    force = solver.force
    ell = np.arcsin(.2)/omega
    exact = 2*radius*momentum+force*C/omega*(ell-np.sin(2*omega*ell)/(2*omega))
    initial = solver.poincare(257)
    assert initial['orbital_angular_momentum_gev_cm'][2] == pytest.approx(exact, abs=1e-10)
    solver.advance(20., .1)
    final = solver.poincare(257)
    np.testing.assert_allclose(final['orbital_angular_momentum_gev_cm'],
                              initial['orbital_angular_momentum_gev_cm'], atol=1e-10)
    np.testing.assert_allclose(final['boost_charge_gev_cm'], initial['boost_charge_gev_cm'], atol=1e-10)


def test_step_budget_is_explicit_and_keeps_a_valid_state():
    solver = ConformalString.radial(225., 774.)
    assert solver.advance(10., .1, max_steps=2) == 'STOP_STEP_BUDGET'
    assert solver.t < 10.
    assert solver.snapshot()['chord_closure_error_cm'] < 1e-10
