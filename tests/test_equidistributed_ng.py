"""Gauge-constrained no-artificial-mass NG finite elements, numerical gates."""
import math
import numpy as np
import pytest

from quirksim.equidistributed_ng import EquidistributedNGString
from quirksim.variational_string import StringParameters
from quirksim.radial_ng_reference import ExactRadialNG,RadialNGParameters


def test_gauge_eliminates_unphysical_interior_point_masses():
    p=StringParameters(nodes=5,initial_length_cm=.1,initial_beta=.4)
    solver=EquidistributedNGString(p)
    assert np.all(solver.masses[1:-1]==0)
    assert solver.gauge_mass==0
    y=solver.initial_state()
    J,b,g,gdot=solver.constraints(y)
    assert J.shape==(3,15)
    assert np.max(np.abs(g))<1e-13
    assert np.max(np.abs(gdot))<1e-13
    energy=solver.action_quantities(y)['energy_gev']
    ref=ExactRadialNG(RadialNGParameters(initial_separation_cm=.1))
    assert energy==pytest.approx(ref.energy_total,abs=1e-10,rel=0)


@pytest.mark.parametrize('nodes',[3,5,7])
def test_gauge_constrained_radial_short_time_matches_exact_ng(nodes):
    p=StringParameters(nodes=nodes,initial_length_cm=.1,initial_beta=.4,
                       b_tesla=0,max_step_ps=.1)
    solver=EquidistributedNGString(p)
    diag,solution=solver.solve(duration_ps=10)
    ref=ExactRadialNG(RadialNGParameters(initial_separation_cm=.1)).state(10)
    x,v=solver._unpack(solution.y[:,-1])
    assert x[-1,0]-x[0,0]==pytest.approx(
        ref['unsigned_separation_cm'],abs=3e-6,rel=0)
    assert diag['max_equal_arclength_deviation_cm']<1e-6
    assert abs(diag['energy_drift_gev'])<1e-5
    assert diag['linear_momentum_drift_gev']<1e-5


def test_magnetic_field_does_not_heat_equal_arclength_string_short_time():
    p=StringParameters(nodes=5,initial_length_cm=.1,initial_beta=.4,
                       b_tesla=2,max_step_ps=.1)
    s=EquidistributedNGString(p)
    y=s.boosted_straight_initial_state(beta_cm_x=.25)
    diag,sol=s.solve(y,duration_ps=4.)
    assert abs(diag['energy_drift_gev'])<1e-5
    assert diag['max_equal_arclength_deviation_cm']<1e-6
    assert diag['max_equal_arclength_rate_residual_cm_per_ps']<1e-6
