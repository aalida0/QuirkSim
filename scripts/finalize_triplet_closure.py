"""Freeze a reproducible closure family and retain unresolved spin sectors."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.angular_radiation import reachability
from quirksim.physical_closure import capture_annihilation_interval,centrifugal_action_floor


def calculate(source,stress_source):
    original=json.loads(gzip.decompress(source.read_bytes()));stress=json.loads(gzip.decompress(stress_source.read_bytes()))
    digest=hashlib.sha256(source.read_bytes()).hexdigest()
    if stress['source_sha256']!=digest or stress['archive_sha256']!=original['archive_sha256']:
        raise ValueError('mismatched provenance')
    model=original['runs'][0]['path']['model'];lam=model['lambda_ev']*1e-9
    if (model['mass_gev'],model['lambda_ev'])!=(225.,774.):raise ValueError('unsupported benchmark')
    width=model['singlet_annihilation_width_gev'];out=[]
    for fit in original['fits']:
        if fit['level']!=8:continue
        subset=[r for r in stress['runs'] if r['case']==fit['case'] and r['gate'] in (1.,10.) and r['loss'] in (0.,1.)]
        intervals=[]
        for r in subset:
            if r['status']=='SECULAR_UPDATE_OUTSIDE_ALLOWED_CENTRAL_ORBITS':
                intervals.append({'lower':0.,'upper':1.});continue
            captured=r['capture_time_ps'] is not None
            age=max(0.,stress['horizon_ps']-r['capture_time_ps']) if captured else 0.
            intervals.append(capture_annihilation_interval(float(captured),age,width))
        generous_energy=fit['fitted_energy_gev']+model['captured_binding_energy_gev']
        budget=reachability(fit['fitted_L_over_hbar'],generous_energy,lam)
        action=centrifugal_action_floor(model['mass_gev'],generous_energy,
                                       budget['deterministic_remaining_L_lower_bound'])
        baseline=next(r['path'] for r in original['runs'] if r['case']==fit['case']
            and r['control']=='fitted_recoil' and r['path']['horizon_ps']==1e13
            and r['path']['model']['epsilon_per_full_period']==.1
            and r['path']['model']['angular_loss_per_effective_emission']==1.)
        captured=baseline['capture_time_ps'] is not None
        age=1e13-baseline['capture_time_ps'] if captured else 0.
        short=next(r['path'] for r in original['runs'] if r['case']==fit['case']
            and r['control']=='fitted_recoil' and r['path']['horizon_ps']==1e3
            and r['path']['model']['epsilon_per_full_period']==.1
            and r['path']['model']['angular_loss_per_effective_emission']==1.)
        short_captured=short['capture_time_ps'] is not None
        short_age=1e3-short['capture_time_ps'] if short_captured else 0.
        out.append({'case':fit['case'],'point':fit['point'],'mode':fit['mode'],'seed':fit['seed'],
            'baseline_status':baseline['status'],'capture_time_ps':baseline['capture_time_ps'],
            'one_ns_post_handoff':{'status':short['status'],
                'annihilation_interval':capture_annihilation_interval(float(short_captured),short_age,width)},
            'baseline_annihilation_intervals_by_spin':{str(f):capture_annihilation_interval(float(captured),age,width,f) for f in (None,0.,.25,1.)},
            'closure_family_annihilation_interval':{'lower':min(i['lower'] for i in intervals),'upper':max(i['upper'] for i in intervals)},
            'closure_paths_checked':len(subset),'unit_emission_budget':budget,'hard_annihilation_barrier':action})
    return {'status':'PHENOMENOLOGICAL_CLOSURE_DEFINED_PHYSICAL_INCLUSIVE_RATE_UNDETERMINED',
        'source_sha256':digest,'stress_source_sha256':hashlib.sha256(stress_source.read_bytes()).hexdigest(),
        'archive_sha256':original['archive_sha256'],'cases':out,
        'prescription':{'scope':'post-first-return, both-neutral, vacuum Dirac triplet; central semiclassical reduction',
            'mass_gev':225.,'lambda_tension_ev':774.,'n_ic':2,'alpha_ic':.04,
            'mean_glueball_energy_loss':'epsilon*Lambda per full signed radial period',
            'epsilon_scenarios':[0.,.01,.1,1.],'gate_radius_in_inverse_lambda':[1.,10.],
            'angular_envelopes':['L conserved','maximal downward 1 hbar per effective Lambda emission; not a preferred drift law'],
            'capture':'absorbing matching boundary at E<=8 Lambda and L<=1',
            'capture_threshold_sensitivities':'selected E=1,32 Lambda and L=0,10 cards; see source stress artifact',
            'capture_matching_efficiency':1.,'capture_efficiency_status':'declared idealization, not calculated overlap',
            'capture_spin':'unspecified by default; no automatic 1/4 population',
            'specified_singlet_width_gev':width,'binding_energy_gev':model['captured_binding_energy_gev'],
            'unmodeled_spin_sector':'annihilation probability retained in [0,1] after capture',
            'excited_constant_hazard':0.,'horizon_ps':1e13},
        'release':{'conditional_scenario_tables':True,'unique_physical_annihilation_probability':False,
            'inclusive_population':False,'uncertainty_band_interpretation':False},
        'missing_physical_inputs':['quantum occupation at matching and spin-one annihilation/transition rates',
            'off-gate emission or large-angular/coherent radiation not bounded by this prescription',
            'mixed-charge and pre-first-return evolution not covered by this closure']}


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]/'studies/data'
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=root/'triplet_ensemble_capture_models_2026-11-11.json.gz')
    p.add_argument('--stress',type=Path,default=root/'triplet_capture_stress_2026-11-12.json.gz');p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=calculate(a.source,a.stress);a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print('frozen cases',len(r['cases']))
