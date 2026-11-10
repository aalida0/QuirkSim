"""Physics sanity gates for all weak decay channels and trajectory comparisons."""
import math
from dataclasses import replace
import numpy as np
import pytest
from quirksim.physics import benchmark, ME, MMU
from quirksim.decays import sample_decay, scalar, leptonic_decay_rest
from quirksim.dynamics import Parameters,propagate,sample_plan,HBARC_GEV_CM


@pytest.mark.parametrize('channel',['pion','electron_neutrino','muon_neutrino'])
@pytest.mark.parametrize('mom',[[0,0,0],[100,30,-70]])
def test_all_channels_exact_four_momentum(channel,mom):
    M=225.16;m=225.;rng=np.random.default_rng(42)
    four=sample_decay(channel,mom,M,m,rng)
    total=sum(four.values(),np.zeros(4))
    assert np.allclose(total,[math.hypot(M,np.linalg.norm(mom)),*mom],atol=3e-9,rtol=0)
    for name,v in four.items():
        mass={'daughter':m,'pion':.13957039,'charged_lepton':ME if channel.startswith('electron') else MMU,
              'neutrino':0.}[name]
        assert abs(scalar(v,v)-mass*mass)<1e-7


def test_leptonic_both_endpoints_complete():
    b=benchmark('doublet',150.)
    plan=sample_plan(b['lifetime_ps'],b['widths_gev'],1)
    assert plan[0]['channel']=='electron_neutrino'
    out=propagate(Parameters(representation='doublet',mass_neutral_gev=150.,
                             lambda_ev=2154.,seed=1,duration_ps=600))
    assert out['end_time_ps']==600
    assert len(out['events'])==2
    assert out['charged']==[False,False]
    assert any(e['channel']=='electron_neutrino' for e in out['events'])
    assert all(e['energy_balance_residual_gev'] is not None and
               abs(e['energy_balance_residual_gev'])<2e-8 for e in out['events'])


def test_matched_replay_identical_latent_decays_and_nontrivial_contact():
    base=Parameters(representation='triplet',mass_neutral_gev=225.,lambda_ev=774.,
                    seed=42,duration_ps=1000)
    recoil=propagate(base)
    ctrl=propagate(replace(base,mode='no_recoil'))
    assert recoil['plan']==ctrl['plan']
    assert [x['channel'] for x in recoil['events']]==[x['channel'] for x in ctrl['events']]
    assert recoil['end_time_ps']==ctrl['end_time_ps']==1000
    assert recoil['post_decay_encounter_count']>0
    assert ctrl['post_decay_encounter_count']>0
    assert recoil['closest_post_decay_turning_point_cm']>ctrl['closest_post_decay_turning_point_cm']
    assert ctrl['post_decay_geometric_core_encounters']>recoil['post_decay_geometric_core_encounters']


def test_stable_exact_radial_period_first_return():
    base=Parameters(seed=42,mode='stable',duration_ps=250,max_step_ps=0.5)
    res=propagate(base)
    M=benchmark('triplet',225.)['charged_mass_gev']
    p=M*base.beta/math.sqrt(1-base.beta**2)
    sigma=(base.lambda_ev*1e-9)**2
    period=2*p*6.582119569e-13/sigma
    assert res['turning_points']
    assert abs(res['turning_points'][0]['time_ps']-period)<0.1
    assert res['first_decay_time_ps'] is None
    assert res['closest_post_decay_turning_point_cm'] is None


def test_string_regulator_and_solver_step_convergence():
    p=Parameters(seed=42,mode='recoil',duration_ps=1000)
    a=propagate(p)
    b=propagate(replace(p,max_step_ps=0.5))
    assert abs(a['closest_post_decay_turning_point_cm']-b['closest_post_decay_turning_point_cm'])<5e-9


def test_post_decay_contact_radius_is_not_an_annihilation_probability():
    p=Parameters(lambda_ev=1000,mode='stable',duration_ps=400)
    r=propagate(p)
    assert r['core_radius_cm']==pytest.approx(HBARC_GEV_CM*1e6)
    assert r['post_decay_geometric_core_encounters']==0
    assert 'no radiation or annihilation' in r['caveats']


def test_directional_kinematics_for_muon_decay_positive():
    rng=np.random.default_rng(197)
    for _ in range(20):
        x=leptonic_decay_rest(150.28,150.,MMU,rng)
        assert x['neutrino'][0]>=0
        assert x['charged_lepton'][0]>=MMU
