"""Reuse the existing orbit-average solver with explicit literature model cards.

This mean-loss/first-return eligibility study does not create radiation events
or assign an annihilation probability to an unmodified classical trajectory.
"""
import argparse
import json
import math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.physics import benchmark
from quirksim.radiation_annihilation import OrbitAveragedModel,linear_cooling_time_ps


def eligibility(lower_mm,upper_mm,radius_mm):
    if (not all(map(math.isfinite,[lower_mm,upper_mm,radius_mm]))
            or lower_mm<0 or upper_mm<lower_mm or radius_mm<=0):
        raise ValueError('ordered nonnegative bounds and positive radius required')
    return ([0.,0.] if lower_mm>radius_mm else [1.,1.] if upper_mm<radius_mm else [0.,1.])


def calculate(window_result):
    lam=774.;b=benchmark('triplet',225.);cards=[]
    for sector,mass,charge in [('charged_stable_reference',b['charged_mass_gev'],1.),
                               ('neutral_reference',225.,0.)]:
        energy=2*mass*(1/math.sqrt(1-.4**2)-1)
        for epsilon in (0.,.01,.1,1.):
            for mean_multiplier in (1,2):
                # Existing full-signed-period implementation is reused. The
                # second card doubles the mean emission count per period to
                # expose the convention ambiguity without modifying dynamics.
                model=OrbitAveragedModel(mass,lam,epsilon,mean_multiplier*lam*1e-9,
                    0.,0.,b['lifetime_ps'],
                    'Li et al. 2311.15486v2 Eq7 mean loss; declared full-period convention',
                    charge_magnitude=charge,include_photons=bool(charge))
                cooling=linear_cooling_time_ps(model,energy,8*lam*1e-9)
                rates=model.rates(energy)
                cards.append({'sector':sector,'epsilon':epsilon,
                    'mean_emission_count_multiplier':mean_multiplier,
                    'initial_pair_excess_energy_gev':energy,
                    'stop_energy_gev':8*lam*1e-9,
                    'stop_convention':'declared end of linear reference; not capture or annihilation',
                    'full_signed_period_ps':rates['full_signed_period_ps'],
                    'mean_loss_per_full_signed_period_gev':mean_multiplier*epsilon*lam*1e-9,
                    'reference_cooling_time_ps':cooling if math.isfinite(cooling) else None,
                    'status':'FINITE_REFERENCE_COOLING_TIME' if math.isfinite(cooling) else 'NO_COOLING_IN_THIS_REFERENCE',
                    'comparison_to_charged_endpoint_weak_lifetime':cooling/b['lifetime_ps'] if math.isfinite(cooling) else None,
                    'applicability':'radial orbit with persistent charge and all encounters eligible; not the decaying G4 trajectory'})
    data=json.loads(window_result.read_text());gates=[]
    for c in data['cases']:
        lower=c['speed_cone']['minimum_separation_lower_bound_mm']
        upper=c['speed_cone']['minimum_separation_upper_bound_mm']
        for multiple in (1.,2.,8.):
            radius=multiple*1.973269804e-4/lam
            eligible=eligibility(lower,upper,radius)
            gates.append({'case':c['name'],'radius_multiple':multiple,'radius_mm':radius,
                          'frame':'lab-frame saved-history gate diagnostic; source CM criterion not yet matched',
                          'radius_choice':'source-motivated scale' if multiple==1 else 'additional sensitivity choice',
                          'eligibility_bounds':eligible,
                          'epsilon_scenarios':[{'epsilon':e,'probability_bounds_for_one_assumed_trial':[e*x for x in eligible]}
                                               for e in (0.,.01,.1,1.)],
                          'meaning':'GEOMETRIC_ELIGIBILITY_ONLY; ONE_TRIAL_IS_A_DECLARED_MODEL_CHOICE; NO_ANNIHILATION_PROBABILITY'})
    return {'scope':'LITERATURE_MOTIVATED_REFERENCE_MODEL_AND_FIRST_RETURN_GATE',
            'sources':{'mean_loss':'2311.15486v2 Eq7',
                       'full_signed_period':'2504.02940v2 Eq6, T=4p/sigma',
                       'zero_glueball_baseline':'2504.02940v2 Sec2.1',
                       'source_scale_warning':'2504.02940 uses sigma=3.6 Lambda_prime^2; our Lambda is sqrt(sigma)'},
            'charged_endpoint_weak_lifetime_ps':b['lifetime_ps'],'cards':cards,'first_return_gates':gates,
            'limitations':['mean emission energy is a loss parameter, not an on-shell glueball spectrum',
                           'no recoil or radiation feedback on trajectories',
                           'no transfer of scalar annihilation or persistent-charge cooling to neutral daughters',
                           'period and emission-count choices are explicit scenarios, not statistical errors',
                           'spin relaxation, angular-momentum loss and capture remain model choices'],
            'physical_annihilation_probability':None}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('window_result',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=calculate(a.window_result);a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    finite=[c['reference_cooling_time_ps'] for c in r['cards'] if c['reference_cooling_time_ps'] is not None]
    print('cards',len(r['cards']),'minimum radial reference cooling time ps',min(finite))
