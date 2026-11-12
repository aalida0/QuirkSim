"""Bridge-level tests for real 4-momentum initialization and topology."""
import math
import numpy as np
import pytest
from quirksim.production_dynamics import (
    ProductionParameters, validated_endpoints, simulate_produced
)
from quirksim.physics import benchmark


def record(pdg1=9000111,pdg2=-9000112):
    cm=benchmark('triplet',225.)['charged_mass_gev']
    data=[]
    for pid,pvec in ((pdg1,[100.,10.,230.]),(pdg2,[-80.,-5.,-160.])):
        mass=225. if abs(pid)==9000112 else cm
        energy=math.hypot(mass,math.sqrt(sum(z*z for z in pvec)))
        q=(1 if abs(pid)==9000111 else -1 if abs(pid)==9000113 else 0)
        if pid<0:q=-q
        data.append({'pdg_id':pid,'charge_e':q,'quirk_number':1 if pid>0 else -1,
                     'mass_gev':mass,'p4_gev':[energy,*pvec]})
    return {'schema_version':1,'event_index':1,'stage':'PYTHIA8_LO_QCD_ISR_FSR',
            'pair_pt_gev':math.hypot(20.,5.),'endpoints':data}


def test_charged_neutral_is_not_two_charged():
    ev=record()
    validated_endpoints(ev)
    out=simulate_produced(ev,ProductionParameters(mode='stable',duration_ps=.3,
                                                  max_step_ps=.1,b_tesla=0))
    assert out['initial_charged_endpoints']==1
    assert out['final_charged_endpoints']==1
    assert out['n_weak_decays']==0
    assert np.allclose(out['initial_p4_gev'],[x['p4_gev'] for x in ev['endpoints']])
    assert out['pair_pt_gev']==ev['pair_pt_gev']


def test_two_charged_charge_signs():
    ev=record(9000111,-9000111)
    out=simulate_produced(ev,ProductionParameters(mode='stable',duration_ps=.3,
                                                  max_step_ps=.1,b_tesla=2))
    assert out['initial_charges_e']==[1,-1]
    assert out['initial_charged_endpoints']==2


def test_q2_charged_sign():
    ev=record(9000113,-9000112)
    validated_endpoints(ev)
    out=simulate_produced(ev,ProductionParameters(mode='stable',duration_ps=.3,
                                                  max_step_ps=.1,b_tesla=0))
    assert out['initial_charges_e']==[-1,0]


def test_wrong_quirk_number_fails():
    ev=record()
    ev['endpoints'][1]['quirk_number']=1
    with pytest.raises(ValueError, match='quirk and antiquark'):
        validated_endpoints(ev)


def test_wrong_charge_fails():
    ev=record()
    ev['endpoints'][0]['charge_e']=-1
    with pytest.raises(ValueError,match='Electric charge'):
        validated_endpoints(ev)


def test_offshell_fails():
    ev=record()
    ev['endpoints'][0]['p4_gev'][0]+=50
    with pytest.raises(ValueError,match='on shell'):
        validated_endpoints(ev)


def test_incompatible_weak_mass_scheme_fails():
    ev=record()
    ev['endpoints'][0]['mass_gev']=226.
    with pytest.raises(ValueError,match='mass'):
        validated_endpoints(ev)
