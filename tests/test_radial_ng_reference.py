"""Independent exact radial (1+1D) massive-endpoint Nambu--Goto gates."""
import math
import numpy as np
import pytest

from quirksim.radial_ng_reference import RadialNGParameters,ExactRadialNG
from quirksim.dynamics import C_CM_PER_PS


def test_radial_first_crossing_and_period_over_thousand_ps():
    p=RadialNGParameters()
    s=ExactRadialNG(p)
    assert 190<s.first_crossing_ps<260
    assert 190<s.crossing_interval_ps<260
    assert s.state(0)['unsigned_separation_cm']==pytest.approx(
        p.initial_separation_cm,rel=1e-11,abs=1e-12)
    assert s.state(s.first_crossing_ps)['unsigned_separation_cm']<1e-10
    assert s.state(s.first_crossing_ps-.01)['crossings']==0
    assert s.state(s.first_crossing_ps+.01)['crossings']==1
    assert s.state(s.first_crossing_ps+s.crossing_interval_ps+.01)['crossings']==2
    assert s.state(1000)['crossings']>=4
    assert s.state(1000)['time_ps']==1000


def test_radial_noether_quantities_at_many_crossings():
    s=ExactRadialNG()
    times=list(np.linspace(0,1000,101))
    times += [s.first_crossing_ps+k*s.crossing_interval_ps
              for k in range(12)]
    values=[s.state(t) for t in times]
    assert max(abs(e['total_energy_gev']-s.energy_total)
               for e in values)<1e-10
    assert all(e['total_linear_momentum_gev']==0 for e in values)
    assert all(e['total_angular_momentum_hbar']==0 for e in values)
    assert all(e['unsigned_separation_cm']>=0 for e in values)
    assert all(e['proper_time_each_ps']>=0 for e in values)


def test_exact_radial_velocity_force_and_proper_time_derivatives():
    s=ExactRadialNG()
    for t in (11.,31.,112.,300.,430.,775.):
        if min(abs(t-s.first_crossing_ps-k*s.crossing_interval_ps)
               for k in range(5))<.05:
            continue
        a=s.state(t-1e-4);b=s.state(t);c=s.state(t+1e-4)
        numerical_beta=(c['right_position_cm']-a['right_position_cm'])/(2e-4*C_CM_PER_PS)
        numerical_force=(c['right_momentum_gev']-a['right_momentum_gev'])/(2e-4)
        numerical_proper=(c['proper_time_each_ps']-a['proper_time_each_ps'])/2e-4
        assert numerical_beta==pytest.approx(b['right_beta'],rel=4e-6,abs=2e-6)
        assert numerical_force==pytest.approx(
            -s.force if b['signed_relative_position_cm']>0 else s.force,
            rel=2e-6)
        assert numerical_proper==pytest.approx(
            s.mass/b['endpoint_energy_each_gev'],rel=2e-6)


def test_radial_phase_invariant_from_gauge_mesh():
    """This exact reduced solution has no interior gauge or artificial masses."""
    s=ExactRadialNG(RadialNGParameters(initial_separation_cm=.1))
    reference=s.state(1000)
    for sampling in (11,31,101,1001):
        out=s.summary(1000,sampling)
        assert out['crossings_at_end']==reference['crossings']
        assert out['final_unsigned_separation_cm']==pytest.approx(
            reference['unsigned_separation_cm'])
        assert out['max_relative_energy_residual']<1e-12


def test_invalid_radial_reference_inputs_rejected():
    with pytest.raises(ValueError):
        RadialNGParameters(lambda_ev=0)
    with pytest.raises(ValueError):
        RadialNGParameters(initial_beta=1)
    with pytest.raises(ValueError):
        ExactRadialNG(RadialNGParameters(initial_beta=0,initial_separation_cm=0))
    with pytest.raises(ValueError):
        ExactRadialNG().state(-.1)


@pytest.mark.parametrize('t',[0.,10.,123.,217.96061652,220.03789114,250.,701.,1000.])
def test_radial_gauge_reparametrization_preserves_physical_charges(t):
    """Irregular interior labellings cannot change 1D string observables."""
    s=ExactRadialNG()
    gauges=[np.linspace(0.,1.,7),
            np.linspace(0.,1.,7)**2,
            np.sqrt(np.linspace(0.,1.,7)),
            np.linspace(0.,1.,13)**3]
    reference=s.state(t)
    for f in gauges:
        sl=s.gauge_slice(t,f)
        assert sl['total_energy_gev']==pytest.approx(
            reference['total_energy_gev'],abs=2e-10,rel=0)
        assert sl['string_energy_gev']==pytest.approx(
            reference['string_energy_gev'],abs=2e-10,rel=0)
        assert sl['endpoint_momentum_gev']==0
        assert sl['string_longitudinal_momentum_gev']==0
        assert sl['vertex_positions_cm'][0]==pytest.approx(
            reference['left_position_cm'])
        assert sl['vertex_positions_cm'][-1]==pytest.approx(
            reference['right_position_cm'])
        assert sl['vertex_betas'][0]==pytest.approx(-reference['right_beta'])
        assert sl['vertex_betas'][-1]==pytest.approx(reference['right_beta'])
        assert max(map(abs,sl['vertex_betas']))<1


def test_relabelling_can_cure_radial_mesh_condition_number_without_changing_energy():
    s=ExactRadialNG()
    time=217.9606165202423
    u_bad=np.array([0.,.01,.5,.75,1.])
    u_good=np.linspace(0,1,len(u_bad))
    bad=s.gauge_slice(time,u_bad)
    good=s.gauge_slice(time,u_good)
    assert bad['max_element_ratio']>20
    assert good['max_element_ratio']==pytest.approx(1.)
    assert good['total_energy_gev']==pytest.approx(
        bad['total_energy_gev'],abs=1e-10,rel=0)
    assert s.state(time)['unsigned_separation_cm']>.05


def test_relabel_rejects_invalid_monotone_gauges():
    s=ExactRadialNG()
    with pytest.raises(ValueError):
        s.gauge_slice(10,[0,0,1])
    with pytest.raises(ValueError):
        s.gauge_slice(10,[0,.3,.2,1])
    with pytest.raises(ValueError):
        s.gauge_slice(10,[0,.5,.9])
