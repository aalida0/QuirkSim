"""Short-time equivalence checks for bounded-speed, massless-gauge NG string."""
import numpy as np
import pytest
from quirksim.equidistributed_ng import EquidistributedNGString
from quirksim.variational_string import StringParameters
from quirksim.variational_rapidity import physical_to_rapidity,rapidity_to_physical
from scripts.rapidity_constrained_contact import run_case


def test_rapidity_constrained_coordinate_roundtrip_for_zero_mass_gauge():
    s=EquidistributedNGString(StringParameters(nodes=5,
                           initial_beta=.4,initial_length_cm=.1))
    state=s.boosted_straight_initial_state(beta_cm_x=.25)
    w=physical_to_rapidity(state,s.n)
    recovered=rapidity_to_physical(w,s.n)
    assert np.max(np.abs(recovered-state))<1e-12


def test_short_field_trajectory_matches_old_coordinate_physical_observables():
    # Verify with the exact parameters used in the standalone cross-check.
    p=StringParameters(nodes=5,initial_beta=.4,initial_length_cm=.1,
                       b_tesla=2.,regularizer_cm=1e-9,max_step_ps=.1)
    s=EquidistributedNGString(p)
    y=s.boosted_straight_initial_state(beta_cm_x=.25)
    base,sol=s.solve(y,duration_ps=4.)
    x,v=s._unpack(sol.y[:,-1])
    rapid=run_case(nodes=5,duration_ps=4,regulator_cm=1e-9,
                   max_step_ps=.1,pinch_fraction=4e-4)
    assert rapid['status']=='REACHED_HORIZON'
    assert rapid['max_vertex_beta']<1
    assert rapid['max_gauge_length_deviation_cm']<1e-7
    assert rapid['energy_drift_gev']==pytest.approx(base['energy_drift_gev'],abs=1e-7)
    assert abs(rapid['endpoint_chord_cm']-
              float(np.linalg.norm(x[-1]-x[0])))<1e-7


def test_solver_returns_explicit_unvalidated_physics_scope():
    result=run_case(nodes=3,duration_ps=.05,regulator_cm=1e-9,
                    pinch_fraction=4e-4,max_step_ps=.05)
    assert result['status']=='REACHED_HORIZON'
    assert 'NOT validated continuum NG' in result['physics_scope']
