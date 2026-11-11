"""Rapidity reparameterization equivalence and enforced causal vertex speeds."""
import numpy as np
import pytest
from quirksim.variational_string import StringParameters, VariationalString
from quirksim.variational_rapidity import (physical_to_rapidity,
    rapidity_to_physical,rhs_rapidity,solve_rapidity)


def test_rapidity_conversion_exact_near_relativistic_speed():
    s=VariationalString(StringParameters(nodes=5,initial_beta=.8))
    y=s.initial_state(transverse_velocity=.02,middle_offset_cm=.001)
    w=physical_to_rapidity(y,s.n)
    z=rapidity_to_physical(w,s.n)
    assert np.allclose(z,y,atol=1e-15,rtol=0)
    w[3*s.n]=1e3
    z=rapidity_to_physical(w,s.n)
    _,v=s._unpack(z)
    assert np.max(np.sum(v*v,axis=1))<1


def test_rapidity_physics_matches_velocity_short_time():
    p=StringParameters(nodes=5,initial_beta=.3,
                       initial_length_cm=.1,duration_ps=.2,max_step_ps=.025)
    s=VariationalString(p)
    y=s.initial_state(middle_offset_cm=5e-5,transverse_velocity=.02)
    old,sol0=s.solve(y)
    fresh,sol1=solve_rapidity(s,y,.2,max_step_ps=.025)
    assert abs(fresh['energy_drift_gev'])<1e-6
    assert fresh['momentum_drift_gev']<1e-6
    assert np.max(np.abs(sol0.y[:,-1]-sol1.y[:,-1]))<1e-7
    assert fresh['max_vertex_beta']<1
