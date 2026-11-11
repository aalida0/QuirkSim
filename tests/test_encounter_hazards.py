"""Hazard envelope is an accounting layer, not a physics prediction."""
import math
import pytest

from quirksim.encounter_hazards import HazardModel,evaluate,scenario_envelope


def test_survival_envelope_and_pre_decay_conditional_fraction():
    events=[{'time_ps':1.,'separation_cm':1e-9},
            {'time_ps':2.,'separation_cm':1e-9},
            {'time_ps':3.,'separation_cm':1e-9}]
    model=HazardModel(774,225,annihilation_per_geometric_contact=.5,
                      glueball_emission_per_contact=.5)
    r=evaluate(events,2.5,model)
    assert r['contacts']==3
    assert r['pre_weak_decay_contacts']==2
    assert r['conditional_survival_to_first_weak_decay']==pytest.approx(.25)
    assert r['conditional_survival_to_observation_end']==pytest.approx(.125)
    assert r['expected_glueball_radiation_gev_without_backreaction']==pytest.approx(
        (1+.5+.25)*.5*774e-9)
    assert r['status']=='ASSUMPTION_ENVELOPE_ONLY'


def test_annihilation_brackets_cover_zero_to_certain():
    events=[{'time_ps':1.,'separation_cm':0.}]
    none=evaluate(events,None,HazardModel(774,225,0.,0.))
    certain=evaluate(events,None,HazardModel(774,225,1.,1.))
    assert none['conditional_survival_to_observation_end']==1
    assert certain['conditional_survival_to_observation_end']==0
    assert none['expected_glueball_radiation_gev_without_backreaction']==0
    assert certain['expected_glueball_radiation_gev_without_backreaction']==pytest.approx(774e-9)


def test_infrared_and_short_distance_scales_differ_by_many_orders():
    l=HazardModel(774,225,radius_definition='one_over_lambda')
    m=HazardModel(774,225,radius_definition='one_over_mass')
    assert l.radius_cm/m.radius_cm>1e8
    assert evaluate([{'time_ps':1,'separation_cm':1e-9}],None,l)['contacts']==1
    assert evaluate([{'time_ps':1,'separation_cm':1e-9}],None,m)['contacts']==0


def test_scenario_envelope_strict_validation():
    r=scenario_envelope([],None,774,225)
    assert len(r)==3*4*3
    assert all(x['conditional_survival_to_observation_end']==1 for x in r)
    with pytest.raises(ValueError):
        HazardModel(774,225,annihilation_per_geometric_contact=1.1)
    with pytest.raises(ValueError):
        evaluate([{'time_ps':1.,'separation_cm':-1.}],None,HazardModel(774,225))
