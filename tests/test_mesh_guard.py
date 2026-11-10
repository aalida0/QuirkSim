"""Strict worldsheet mesh safety prevents interpreting extrapolated trajectories."""
import numpy as np
import pytest

from quirksim.mesh_guard import (
    MeshQualityViolation,make_mesh_stop_events,mesh_metrics,
    require_good_mesh,is_mesh_acceptable)
from quirksim.variational_string import StringParameters,VariationalString
from quirksim.variational_rapidity import solve_rapidity


def test_quality_events_trigger_on_synthetic_bad_elements():
    nodes=5
    reg=1e-9
    state=np.zeros(6*nodes+2)
    state[:3*nodes]=np.array([[0.,0,0],[.025,0,0],
                               [.05,0,0],[.075,0,0],[.1,0,0]]).ravel()
    require_good_mesh(state,nodes,reg)
    assert is_mesh_acceptable(state,nodes,reg)
    events=make_mesh_stop_events(nodes,reg)
    assert events[0](0.,state)>0
    assert events[1](0.,state)>0
    bad=state.copy()
    bad[3]=1e-11
    assert not is_mesh_acceptable(bad,nodes,reg)
    assert events[0](0.,bad)<0
    assert events[1](0.,bad)<0
    with pytest.raises(MeshQualityViolation):
        require_good_mesh(bad,nodes,reg)


def test_strict_rapidity_rejects_malformed_radial_mesh_before_first_crossing():
    p=StringParameters(nodes=5,initial_beta=.4,initial_length_cm=.1,
                        b_tesla=0,max_step_ps=.5)
    s=VariationalString(p)
    with pytest.raises(MeshQualityViolation,match='element_imbalance'):
        solve_rapidity(s,s.initial_state(),duration_ps=225,
                       max_step_ps=.5,strict_mesh=True)
