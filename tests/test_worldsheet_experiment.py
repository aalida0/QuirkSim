"""Event-driven weak decay pilot smoke gates; does not certify physical worldsheet."""
import math
from quirksim.variational_string import StringParameters
from quirksim.worldsheet_experiment import simulate


def test_early_stochastic_triplet_decay_has_canonical_conservation():
    # Seed 25 happens to select a rare short lifetime for one endpoint.
    p=StringParameters(nodes=5,b_tesla=2,initial_beta=.4,
                       initial_length_cm=.1,max_step_ps=.05)
    r=simulate(seed=25,params=p,duration_ps=2.,beta_cm_x=.25)
    assert r['elapsed_ps']<=2.
    assert r['radiation_annihilation_hazard_envelope_summary']['scenarios']==36
    assert r['status'] in ('FINITE_ELEMENT_PROTOTYPE_ONLY','BLOCKED_PHYSICS_GATE')
    if r['status']=='FINITE_ELEMENT_PROTOTYPE_ONLY':
        assert len(r['decays'])>=1
        assert not r['charged_endpoints'][r['decays'][0]['endpoint_label']]
        assert abs(r['decays'][0]['energy_residual_gev'])<1e-4
        assert max(abs(e) for e in r['string_energy_drift_by_free_segment_gev'])<1e-3
    else:
        assert r['failure_reason'] is not None


def test_short_stable_interval_does_not_make_up_an_annihilation_rate():
    p=StringParameters(nodes=3,b_tesla=0,max_step_ps=.05)
    r=simulate(seed=10,params=p,duration_ps=.05,beta_cm_x=0)
    assert not r['decays']
    assert r['post_decay_closest_cm'] is None
    assert r['radiation_annihilation_hazard_envelope_summary']['possible_survival_to_weak_range']==[None,None]
