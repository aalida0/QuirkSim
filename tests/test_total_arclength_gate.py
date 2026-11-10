"""Strict common TOTAL arclength event does not depend on mesh count or clustering."""
import numpy as np
import pytest
from quirksim.mesh_guard import (
    total_string_arclength_cm,make_total_arclength_stop_event,
    make_mesh_stop_events,MeshQualityViolation)
from quirksim.equidistributed_ng import EquidistributedNGString
from quirksim.variational_string import StringParameters


@pytest.mark.parametrize('n',[5,7,11])
@pytest.mark.parametrize('exponent',[1.,1.3,1.8])
def test_total_arclength_pinch_event_mesh_and_gauge_invariant(n,exponent):
    cutoff=4e-5
    state=np.zeros(6*n+2)
    u=np.linspace(0,1,n)**exponent
    state[:3*n]=np.column_stack([cutoff*u,np.zeros(n),np.zeros(n)]).ravel()
    event=make_total_arclength_stop_event(n,cutoff)
    assert total_string_arclength_cm(state,n)==pytest.approx(cutoff,abs=2e-20)
    assert event(0,state)==pytest.approx(0.,abs=2e-20)
    assert event.terminal
    assert event.direction==-1
    larger=state.copy()
    larger[:3*n]*=1.001
    smaller=state.copy()
    smaller[:3*n]*=.999
    assert event(0,larger)>0
    assert event(0,smaller)<0


def test_longtime_constrained_solver_fail_closed_on_total_arclength():
    params=StringParameters(nodes=5,initial_length_cm=.1,
                            initial_beta=.4,b_tesla=2,
                            max_step_ps=.5,regularizer_cm=1e-9)
    solver=EquidistributedNGString(params)
    state=solver.boosted_straight_initial_state(beta_cm_x=.25)
    with pytest.raises(MeshQualityViolation,match='common total-arclength cutoff'):
        solver.solve(state,duration_ps=480,strict_mesh=True,
                     pinch_fraction_of_initial_length=4e-4)


@pytest.mark.parametrize('invalid',[0.,-1.,float('nan'),float('inf')])
def test_invalid_total_arclength_cuts_rejected(invalid):
    with pytest.raises(ValueError):
        make_total_arclength_stop_event(5,invalid)
