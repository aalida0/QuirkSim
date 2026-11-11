"""Gauss Nambu-Goto action and canonical Euler-Lagrange differentiation."""
import numpy as np
import pytest
from quirksim.variational_string import StringParameters
from quirksim.equidistributed_ng import EquidistributedNGString
from quirksim.ng_quadrature import GaussQuadratureNGString
from quirksim.dynamics import C_CM_PER_PS


def action_L(s,y):
    q=s.action_quantities(y)
    _,v=s._unpack(y)
    return float(np.sum(q['canonical_vertex_momenta_gev']*v)
                 -q['energy_gev'])


def test_order_one_gauss_matches_midpoint_semidiscrete_action():
    p=StringParameters(nodes=5,initial_length_cm=.1,initial_beta=.2)
    baseline=EquidistributedNGString(p)
    gauss=GaussQuadratureNGString(p,1)
    y=baseline.initial_state()
    x,v=baseline._unpack(y)
    v[1]=[.04,.08,.01]
    v[2]=[.0,.02,.02]
    v[3]=[-.04,-.02,-.01]
    x[2,1]=.001
    b=baseline.action_quantities(y,True)
    g=gauss.action_quantities(y,True)
    for k in ('energy_gev','momentum_gev','canonical_vertex_momenta_gev',
              'mass_matrix','forces_gev_per_ps','convective_gev_per_ps'):
        np.testing.assert_allclose(g[k],b[k],rtol=1e-9,atol=1e-8)


@pytest.mark.parametrize('order',[2,3])
def test_gauss_canonical_momenta_are_action_velocity_derivatives(order):
    p=StringParameters(nodes=5,initial_beta=.2)
    s=GaussQuadratureNGString(p,order)
    y=s.initial_state()
    x,v=s._unpack(y)
    x[2,1]=.0003
    v[1]=[.01,.04,-.02]
    v[2]=[-.01,-.05,.01]
    v[3]=[-.02,.02,.015]
    result=s.action_quantities(y,True)
    eps=1e-6
    for k in (0,2,4,6,8,11,14):
        yp=y.copy();ym=y.copy()
        yp[3*s.n+k]+=eps;ym[3*s.n+k]-=eps
        measured=(action_L(s,yp)-action_L(s,ym))/(2*eps)
        assert measured==pytest.approx(
            result['canonical_vertex_momenta_gev'].ravel()[k],abs=6e-8)


@pytest.mark.parametrize('order',[2,3])
def test_gauss_position_forces_follow_action_derivative(order):
    p=StringParameters(nodes=5,initial_beta=.2)
    s=GaussQuadratureNGString(p,order)
    y=s.initial_state()
    x,v=s._unpack(y)
    x[2,1]=.0003
    v[1]=[.01,.04,-.02]
    v[2]=[-.01,-.05,.01]
    v[3]=[-.02,.02,.015]
    result=s.action_quantities(y,True)
    eps=1e-7
    for k in (0,1,3,4,6,7,10,12,14):
        yp=y.copy();ym=y.copy()
        yp[k]+=eps;ym[k]-=eps
        measured=(action_L(s,yp)-action_L(s,ym))/(2*eps)
        assert measured==pytest.approx(
            result['forces_gev_per_ps'].ravel()[k]/C_CM_PER_PS,
            abs=3e-5)


def test_gauss2_string_conserves_energy_for_short_field_run():
    p=StringParameters(nodes=5,initial_beta=.4,initial_length_cm=.1,
                       b_tesla=2,max_step_ps=.1)
    s=GaussQuadratureNGString(p,2)
    y=s.boosted_straight_initial_state(.25)
    summary,solution=s.solve(y,duration_ps=4)
    assert abs(summary['energy_drift_gev'])<1e-6
    assert summary['max_equal_arclength_deviation_cm']<1e-6
    assert summary['max_equal_arclength_rate_residual_cm_per_ps']<1e-6
