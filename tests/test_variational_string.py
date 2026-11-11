"""Action-level and dynamical validation of the spatially discretized string.

Scientific limitation: tested finite-element NG midpoint ACTION, not the
unregulated continuum worldsheet. The discretization and gauge regulator
must converge before weak-decay physics conclusions.
"""
from dataclasses import replace
import math
import numpy as np
import pytest

from quirksim.variational_string import StringParameters,VariationalString
from quirksim.dynamics import C_CM_PER_PS


def lagrangian(s, state):
    q=s.action_quantities(state)
    _,v=s._unpack(state)
    return float(np.sum(q['canonical_vertex_momenta_gev']*v)-q['energy_gev'])


def test_discrete_action_matches_canonical_momenta():
    p=StringParameters(nodes=4,initial_beta=.3)
    s=VariationalString(p)
    y=s.initial_state(transverse_velocity=.08,middle_offset_cm=.0002)
    x,v=s._unpack(y)
    v[1]=[.02,.1,-.03]
    v[2]=[-.02,.04,.05]
    q=s.action_quantities(y)
    delta=1e-6
    for j in (0,3,5,7,11):
        plus=y.copy();minus=y.copy()
        k=3*s.n+j
        plus[k]+=delta;minus[k]-=delta
        numeric=(lagrangian(s,plus)-lagrangian(s,minus))/(2*delta)
        assert numeric==pytest.approx(q['canonical_vertex_momenta_gev'].ravel()[j],
                                      abs=3e-8)


def test_discrete_action_matches_position_force():
    p=StringParameters(nodes=4,initial_beta=.3)
    s=VariationalString(p)
    y=s.initial_state(transverse_velocity=.08,middle_offset_cm=.0002)
    _,v=s._unpack(y)
    v[1]=[.02,.1,-.03]
    v[2]=[-.02,.04,.05]
    q=s.action_quantities(y,derivatives=True)
    delta=1e-7
    for j in (0,4,6,8,11):
        plus=y.copy();minus=y.copy()
        plus[j]+=delta;minus[j]-=delta
        numeric=(lagrangian(s,plus)-lagrangian(s,minus))/(2*delta)
        force=q['forces_gev_per_ps'].ravel()[j]/C_CM_PER_PS
        assert numeric==pytest.approx(force,abs=3e-5)


def test_radial_undecayed_motion_energy_and_momentum():
    p=StringParameters(nodes=5,initial_beta=.4,
                       initial_length_cm=.1,duration_ps=1.)
    s=VariationalString(p)
    result,_=s.solve()
    assert abs(result['energy_drift_gev'])<1e-6
    assert np.linalg.norm(np.array(result['final_momentum_gev'])-
                          result['initial_momentum_gev'])<1e-6
    assert np.linalg.norm(result['final_angular_momentum_hbar'])<1e-2


def test_bent_string_carries_internal_momentum_and_angular_momentum():
    p=StringParameters(nodes=5,initial_beta=.0,duration_ps=.4,
                       initial_length_cm=.1,gauge_mass_fraction=.01)
    s=VariationalString(p)
    y=s.initial_state(middle_offset_cm=.0002)
    _,v=s._unpack(y)
    v[2,2]=.1
    q=s.action_quantities(y)
    assert np.linalg.norm(q['canonical_vertex_momenta_gev'][2])>0
    assert np.linalg.norm(q['angular_momentum_over_hbar'])>0
    result,_=s.solve(y)
    assert abs(result['energy_drift_gev'])<1e-5
    assert np.linalg.norm(np.array(result['final_momentum_gev'])-
                          result['initial_momentum_gev'])<1e-5
    assert np.linalg.norm(np.array(result['final_angular_momentum_hbar'])-
                          result['initial_angular_momentum_hbar'])<2e6


def test_magnetic_force_does_no_work():
    p=StringParameters(nodes=5,initial_beta=.1,
                       initial_length_cm=.1,b_tesla=2,
                       duration_ps=.5)
    s=VariationalString(p)
    y=s.initial_state(transverse_velocity=.03)
    x,v=s._unpack(y)
    q=s.action_quantities(y,derivatives=True)
    # B field should do no instantaneous work on charged endpoints.
    sw=VariationalString(replace(p,b_tesla=0))
    z=sw.action_quantities(y,derivatives=True)
    magnetic=q['forces_gev_per_ps']-z['forces_gev_per_ps']
    assert np.sum(magnetic*v)==pytest.approx(0.,abs=1e-12)


def test_regulator_is_explicit_and_no_weak_decay_claim():
    s=VariationalString(StringParameters())
    r,_=s.solve(duration_ps=.1)
    assert s.gauge_mass>0
    assert any('NOT implemented' in x for x in r['limitations'])
    assert r['model']=='variational_midpoint_Nambu_Goto_finite_elements'
