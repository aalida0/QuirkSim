"""Published straight-string endpoint force and solenoid angular-momentum tests."""
from dataclasses import replace
import math
import numpy as np
import pytest
from quirksim.field_dynamics import (
    FieldParameters, MAGNETIC_FORCE_GEV_PER_PS_T, ng_force,
    initial_state, simulate, string_force_pair)
from quirksim.physics import benchmark, HBAR_GEV_PS


def test_ng_string_radial_exact_equivalence():
    direction=np.array([1.,0.,0.])
    strength=2.345
    for beta in [-.72,0.,.61]:
        force=ng_force(np.array([beta,0.,0.]),direction,strength)
        assert np.allclose(force,[-strength,0,0],atol=1e-14)


def test_ng_string_nonradial_transverse_term():
    f=ng_force(np.array([.4,.3,0]),np.array([1.,0.,0.]),1.)
    assert f[0]==pytest.approx(-math.sqrt(1-.3*.3),rel=1e-14)
    assert f[1]==pytest.approx(-.4*.3/math.sqrt(1-.3*.3),rel=1e-14)


def test_lorentz_force_conversion():
    # 1 GeV charge moving transverse to B=1T curves at p[GeV]/(0.2998 B) meters.
    assert MAGNETIC_FORCE_GEV_PER_PS_T==pytest.approx(
        0.299792458*0.000299792458,rel=1e-14)


def test_coincident_production_no_spurious_initial_orbital_L():
    p=FieldParameters(beta_cm_x=.25)
    mass=benchmark('triplet',225)['charged_mass_gev']
    y=initial_state(p,mass)
    assert np.cross(y[:3]-y[3:6],y[6:9]-y[9:12]).tolist()==[0.,0.,0.]
    assert y[6]>0 and y[9]<0
    assert y[6]+y[9]>0 # boosted pair


def test_magnetic_field_generates_pre_decay_L_from_zero():
    p=FieldParameters(mode='stable',duration_ps=400,max_step_ps=1.0)
    no=simulate(replace(p,b_tesla=0))
    plus=simulate(replace(p,b_tesla=2))
    minus=simulate(replace(p,b_tesla=-2))
    assert no['initial_relative_L_over_hbar']==0
    assert plus['initial_relative_L_over_hbar']==0
    assert no['max_pre_first_decay_L_over_hbar']<1e-3
    assert plus['max_pre_first_decay_L_over_hbar']>1e8
    assert minus['max_pre_first_decay_L_over_hbar']==pytest.approx(
        plus['max_pre_first_decay_L_over_hbar'],rel=1e-8)
    assert plus['turning_points']


def test_symmetric_no_boost_reduces_field_sensitivity():
    p=FieldParameters(mode='stable',duration_ps=300,b_tesla=2, beta_cm_x=0)
    symmetric=simulate(p)
    boosted=simulate(replace(p,beta_cm_x=.25))
    assert boosted['max_pre_first_decay_L_over_hbar']>symmetric['max_pre_first_decay_L_over_hbar']*100


def test_radial_analytic_return_and_ng_vs_central():
    p=FieldParameters(mode='stable',duration_ps=250,b_tesla=0,beta_cm_x=0)
    central=simulate(replace(p,force_model='central_proxy'))
    ng=simulate(replace(p,force_model='ng_endpoint'))
    charged_mass=benchmark('triplet',225)['charged_mass_gev']
    momentum=charged_mass*p.beta_rel/math.sqrt(1-p.beta_rel**2)
    expected=2*momentum*HBAR_GEV_PS/(p.lambda_ev*1e-9)**2
    for sim in (central,ng):
        assert abs(sim['turning_points'][0]['time_ps']-expected)<.1
        energies=sim['energy_proxy_segment_endpoints']
        assert max(energies)-min(energies)<1e-4
    assert abs(central['turning_points'][0]['time_ps']-ng['turning_points'][0]['time_ps'])<1e-6


def test_magnetic_recoil_weak_decay_and_proper_time_matching():
    p=FieldParameters(mode='recoil',duration_ps=1000,seed=42,b_tesla=2)
    physical=simulate(p)
    zero_recoil=simulate(replace(p,mode='no_recoil'))
    assert physical['plan']==zero_recoil['plan']
    assert len(physical['decays'])>=1
    assert len(physical['decays'])==len(zero_recoil['decays'])
    assert [x['channel'] for x in physical['decays']]==[x['channel'] for x in zero_recoil['decays']]
    assert physical['first_decay_L_over_hbar'] is not None
    assert physical['max_pre_first_decay_L_over_hbar']>0
    assert physical['end_time_ps']==1000


def test_ng_endpoint_not_claimed_full_ng_string():
    r=simulate(FieldParameters(mode='stable',duration_ps=5))
    assert any('no worldsheet' in caveat for caveat in r['caveats'])
