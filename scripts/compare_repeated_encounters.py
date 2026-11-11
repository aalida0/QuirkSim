"""Connect uploaded neutral recoil windows to a declared secular orbit model."""
import argparse
import hashlib
import io
import json
import math
from pathlib import Path
import sys
import tarfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.repeated_encounters import EncounterPrescription,evolve,pericenter_gev_inverse
from audit_quirkpheno_return_window import chord_minimum

HBARC_GEV_MM=1.973269804e-13
C_MM_PS=.299792458


def fitted_angular_momentum(case,separation_mm,relative_momentum):
    """Use exact 1D symmetry only for vacuum zero-field CM-collinear controls.

    Dividing tiny interpolation residuals by hbar otherwise produces spurious
    O(100) angular quantum numbers near a symmetry-protected crossing.
    """
    raw=float(np.linalg.norm(np.cross(separation_mm/HBARC_GEV_MM,relative_momentum)))
    correction={'raw_endpoint_L_over_hbar':raw,'symmetry_projection_applied':False}
    if case.get('field_tesla')==0. and case.get('direction_mode')=='collinear_cm':
        theta=case['theta_cm'];axis=np.array([math.sin(theta),0.,math.cos(theta)])
        dx=float(np.linalg.norm(separation_mm-axis*(separation_mm@axis)))
        dp=float(np.linalg.norm(relative_momentum-axis*(relative_momentum@axis)))
        if dx>1e-8 or dp>1e-6:raise ValueError('zero-field collinear control violates expected symmetry')
        correction.update(symmetry_projection_applied=True,transverse_position_residual_mm=dx,
                          transverse_momentum_residual_gev=dp,
                          reason='exact vacuum 1D collinear-control invariant; not a general L rounding threshold')
        return 0.,correction
    return raw,correction


def fit_recoil_windows(archive):
    result=[]
    with tarfile.open(archive) as tar:
        manifest=json.load(tar.extractfile('manifest.json'))
        for case in manifest['cases']:
            if not case['weak_enabled']:continue
            if case.get('window_fit_status','BOTH_NEUTRAL_AT_FIRST_RETURN')!='BOTH_NEUTRAL_AT_FIRST_RETURN':continue
            beta=case.get('pair_boost_z',.25);gamma=1/math.sqrt(1-beta*beta)
            pair=[];momenta=[]
            with np.load(io.BytesIO(tar.extractfile(case['name']+'.npz').read()),allow_pickle=False) as a:
                for q in (0,1):
                    mask=a['q']==q
                    if not np.all(np.abs(a['pdg'][mask])==9900017):raise ValueError('window must be after both weak decays')
                    masses=a['mass'][mask]/1000
                    if np.max(abs(masses-225.))>1e-10:raise ValueError('only neutral 225 GeV benchmark supported')
                    labtime=a['t'][mask]*1000
                    xyz=np.column_stack([a[k][mask] for k in ('x','y','z')])
                    t=gamma*(labtime-beta*xyz[:,2]/C_MM_PS)
                    x=xyz.copy();x[:,2]=gamma*(xyz[:,2]-beta*C_MM_PS*labtime)
                    p=np.column_stack([a[k][mask]/1000 for k in ('px','py','pz')])
                    energy=np.sqrt(masses*masses+np.sum(p*p,axis=1))
                    p[:,2]=gamma*(p[:,2]-beta*energy)
                    if np.any(np.diff(t)<0):raise ValueError('nonmonotone transformed history')
                    keep=np.r_[np.diff(t)>0,True];pair.append((t[keep],x[keep]));momenta.append(p[keep])
            minimum=chord_minimum(pair);time=minimum['time_ps']
            x=[np.array([np.interp(time,pair[i][0],pair[i][1][:,j]) for j in range(3)]) for i in (0,1)]
            p=[np.array([np.interp(time,pair[i][0],momenta[i][:,j]) for j in range(3)]) for i in (0,1)]
            r=(x[0]-x[1])/HBARC_GEV_MM;relative=(p[0]-p[1])/2
            ell,angular_correction=fitted_angular_momentum(case,x[0]-x[1],relative)
            excess=2*(math.hypot(225.,np.linalg.norm(relative))-225.)+(774e-9)**2*np.linalg.norm(r)
            inner=pericenter_gev_inverse(225.,excess,ell,774.)
            result.append({'case':case['name'],'frame':f'fixed production CM beta_z={beta}',
                'mode':case.get('mode'),'point':case.get('point'),'seed':case.get('seed'),'level':case['level'],
                'fit':'instantaneous equal-time endpoint central Hamiltonian reduction; string field energy/angular momentum omitted',
                'time_in_production_cm_ps':time,'separation_um':minimum['minimum_mm']*1000,
                'fitted_energy_gev':excess,'fitted_L_over_hbar':ell,
                'angular_fit_diagnostic':angular_correction,
                'central_inner_turning_point_um':inner*HBARC_GEV_MM*1000,
                'pair_endpoint_residual_momentum_gev':(p[0]+p[1]).tolist()})
    return result


def probabilities(path,spin_fraction,optical_depth,width):
    if path['annihilation_probability_in_declared_model'] is None:return {'status':'TRUNCATED'}
    survive=math.exp(-optical_depth*path['eligible_encounters_in_secular_clock'])
    captured=path['capture_time_ps'] is not None
    age=max(0.,path['horizon_ps']-path['capture_time_ps']) if captured else 0.
    from quirksim.physics import HBAR_GEV_PS
    decay=-math.expm1(-age*width/HBAR_GEV_PS)
    return {'status':'DECLARED_MODEL_ONLY','singlet_fraction_at_capture':spin_fraction,
            'excited_optical_depth_per_eligible_encounter':optical_depth,
            'annihilation_probability':1-survive+(survive*spin_fraction*decay if captured else 0.),
            'capture_probability':survive if captured else 0.}


def paired_control_summary(fits,runs):
    """Within-point paired seed differences; never pool controls as population."""
    metadata={f['case']:f for f in fits};lookup={};groups={}
    for r in runs:
        f=metadata[r['case']]
        if r['control']!='fitted_recoil' or f.get('level')!=8 or f.get('mode') not in ('isotropic','collinear_forward','collinear_backward'):continue
        p=r['path'];m=p['model']
        key=(f['point'],f['seed'],m['epsilon_per_full_period'],m['angular_loss_per_effective_emission'],p['horizon_ps'])
        lookup.setdefault(key,{})[f['mode']]=r['closure_scenarios']
    for key,modes in lookup.items():
        if 'isotropic' not in modes:continue
        point,seed,eps,loss,horizon=key
        for control in ('collinear_forward','collinear_backward'):
            if control not in modes:continue
            for iso,ref in zip(modes['isotropic'],modes[control]):
                if iso['status']!='DECLARED_MODEL_ONLY' or ref['status']!='DECLARED_MODEL_ONLY':continue
                f,h=iso['singlet_fraction_at_capture'],iso['excited_optical_depth_per_eligible_encounter']
                if (f,h)!=(ref['singlet_fraction_at_capture'],ref['excited_optical_depth_per_eligible_encounter']):raise ValueError('mismatched paired closure cards')
                group=(point,control,eps,loss,horizon,f,h)
                groups.setdefault(group,[]).append((seed,iso['annihilation_probability']-ref['annihilation_probability']))
    result=[]
    for key,values in groups.items():
        point,control,eps,loss,horizon,f,h=key;delta=np.array([v for seed,v in values])
        result.append({'point':point,'conditional_direction_control':control,'epsilon':eps,'angular_loss':loss,'horizon_ps':horizon,
            'capture_singlet_fraction':f,'excited_optical_depth':h,'matched_seeds':[seed for seed,v in values],
            'mean_isotropic_minus_control_probability':float(np.mean(delta)),
            'paired_standard_error':float(np.std(delta,ddof=1)/np.sqrt(len(delta))) if len(delta)>1 else None,
            'selection':'both modes neutral-qualified by first return; conditional, not inclusive; unqualified seeds remain in the ensemble manifest'})
    return result


def calculate(archive,compact_result):
    fits=fit_recoil_windows(archive);compact=json.loads(compact_result.read_text())
    with tarfile.open(archive) as tar:input_manifest=json.load(tar.extractfile('manifest.json'))
    expected={'neutral_mass_gev':225.,'n_ic':2,'alpha_ic_hard':.04,'lambda_tension_ev':774.}
    if any(compact['model'].get(k)!=v for k,v in expected.items()):
        raise ValueError('compact-state input does not match the post-weak orbit card')
    width=compact['sum_modeled_two_boson_widths_gev']
    binding=-compact['mesh_checks'][-1]['state']['energy_relative_to_neutral_threshold_gev']
    runs=[]
    for fit in fits:
        for control in ('fitted_recoil','zero_L_action_reference'):
            ell=fit['fitted_L_over_hbar'] if control=='fitted_recoil' else 0.
            for horizon in (1e3,1e13):
                for epsilon in (.01,.1,1.):
                    for loss in (0.,1.):
                        model=EncounterPrescription(epsilon_per_full_period=epsilon,
                            angular_loss_per_effective_emission=loss,captured_binding_energy_gev=binding,
                            singlet_annihilation_width_gev=width,captured_singlet_fraction=1.)
                        path=evolve(fit['fitted_energy_gev'],ell,horizon,model)
                        path.pop('history')
                        scenarios=[probabilities(path,f,h,width) for f in (0.,.25,1.) for h in (0.,.01)]
                        runs.append({'case':fit['case'],'control':control,'path':path,'closure_scenarios':scenarios})
    probe={}
    for loss in (0.,1.):
        r=evolve(.001,100.,1e8,EncounterPrescription(angular_loss_per_effective_emission=loss))
        r.pop('history');probe[str(loss)]=r
    return {'scope':'PHENOMENOLOGICAL_POST_WEAK_REPEATED_ENCOUNTER_MODEL',
        'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
        'compact_result_sha256':hashlib.sha256(compact_result.read_bytes()).hexdigest(),
        'fits':fits,'runs':runs,'synthetic_angular_loss_probe':probe,
        'paired_conditional_direction_comparisons':paired_control_summary(fits,runs),
        'population_warning':'Fits only cover neutral-qualified first-return windows; never infer inclusive survival by dropping other categories',
        'input_population_summary':input_manifest.get('ensemble_summary'),
        'input_conditional_control_summary':input_manifest.get('conditional_control_summary'),
        'limitations':['model begins at a fitted first-return state after both decays',
            'zero-L reference is not a physical no-recoil weak decay control',
            'central reduction omits NG string modes and external fields',
            'radial period is a reference approximation for nonradial orbits',
            'radiation is secular mean energy/action loss, not sampled on-shell particles',
            'directed angular loss, capture and spin reset are declared prescriptions',
            'non-singlet captured sector has no further annihilation model here',
            'pre-handoff radiation and future material/field kicks are absent'],
        'model_independent_annihilation_probability':None}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path)
    p.add_argument('--compact-result',type=Path,default=Path(__file__).resolve().parents[1]/'studies/data/neutral_triplet_quantum_2026-11-11.json')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=calculate(a.archive,a.compact_result);a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    for f in r['fits']:print(f['case'],'fit L/hbar',f['fitted_L_over_hbar'],'pericenter um',f['central_inner_turning_point_um'])
    print('declared orbit scenarios',len(r['runs']))
