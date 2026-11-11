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


def fit_recoil_windows(archive):
    result=[];beta=.25;gamma=1/math.sqrt(1-beta*beta)
    with tarfile.open(archive) as tar:
        manifest=json.load(tar.extractfile('manifest.json'))
        for case in manifest['cases']:
            if not case['weak_enabled']:continue
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
            ell=float(np.linalg.norm(np.cross(r,relative)))
            excess=2*(math.hypot(225.,np.linalg.norm(relative))-225.)+(774e-9)**2*np.linalg.norm(r)
            inner=pericenter_gev_inverse(225.,excess,ell,774.)
            result.append({'case':case['name'],'frame':'fixed production CM beta_z=0.25',
                'fit':'instantaneous equal-time endpoint central Hamiltonian reduction; string field energy/angular momentum omitted',
                'time_in_production_cm_ps':time,'separation_um':minimum['minimum_mm']*1000,
                'fitted_energy_gev':excess,'fitted_L_over_hbar':ell,
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


def calculate(archive,compact_result):
    fits=fit_recoil_windows(archive);compact=json.loads(compact_result.read_text())
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
