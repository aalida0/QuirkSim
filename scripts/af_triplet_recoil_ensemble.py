"""Bounded kinematic-grid pilot with exact conditional pion-direction controls."""
import argparse
import hashlib
import io
import json
import math
from pathlib import Path
import re
import tarfile
import numpy as np
from af_quirkpheno_weak_pilot import prepare as prepare_sampled,execute,M_PARENT,ASSETS
from audit_quirkpheno_return_window import chord_minimum

POINTS=[{'name':'benchmark','beta':.4,'theta':math.pi/2,'field':2.},
        {'name':'tilted','beta':.4,'theta':math.pi/3,'field':2.},
        {'name':'faster','beta':.6,'theta':math.pi/2,'field':2.},
        {'name':'zero_field','beta':.4,'theta':math.pi/2,'field':0.}]


def prepare(out,seeds=8):
    base=prepare_sampled(out,sampled_pion_seeds=seeds)
    templates={level:(out/f'sampled_stable_level{level}.mac').read_text() for level in (7,8)}
    plans={c['seed']:c['plan'] for c in base['cases'] if c['weak_enabled'] and c['mode'].startswith('sampled_seed')}
    for c in base['cases']:
        for suffix in ('.mac','.card','.plan'):
            path=out/(c['name']+suffix)
            if path.exists():path.unlink()
    cases=[]
    for point in POINTS:
        momentum=M_PARENT*point['beta']/math.sqrt(1-point['beta']**2)
        energy=math.hypot(M_PARENT,momentum);boost=.25;gamma=1/math.sqrt(1-boost**2)
        px=momentum*math.sin(point['theta']);pz=momentum*math.cos(point['theta'])
        first_return=2*momentum*6.582119569e-13/(774e-9)**2*gamma
        def add(seed,mode,level):
            name=f"{point['name']}_{mode}_seed{seed}_level{level}"
            card=out/(name+'.card');card.write_text(f"world 100 100 G4_Galactic\nfield {point['field']} 0 100 0 100\n")
            macro=templates[level].replace('sampled_stable_level'+str(level),name)
            macro=re.sub(r'^/qsim/gun/p1 .*$',f'/qsim/gun/p1 {px:.17g} 0 {gamma*(pz+boost*energy):.17g} GeV',macro,flags=re.M)
            macro=re.sub(r'^/qsim/gun/p2 .*$',f'/qsim/gun/p2 {-px:.17g} 0 {gamma*(-pz+boost*energy):.17g} GeV',macro,flags=re.M)
            macro=re.sub(r'^/qsim/maxTime .*$',f'/qsim/maxTime {1.4*first_return/1000:.17g} ns',macro,flags=re.M)
            path=out/(name+'.mac');path.write_text(macro)
            c={'name':name,'mode':mode,'level':level,'seed':seed,'point':point['name'],
               'macro':str(path),'macro_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
               'weak_enabled':mode!='stable','endpoint':-1,'clock_mode':'localized',
               'direction_mode':'collinear_cm' if mode.startswith('collinear') else 'planned',
               'collinear_sign':-1 if mode=='collinear_backward' else 1,
               'field_tesla':point['field'],'beta_each_cm':point['beta'],'theta_cm':point['theta'],
               'pair_boost_z':boost,'nominal_stop_ps':1.4*first_return,'predicted_first_return_ps':first_return,
               'max_boost':.01/2**level,'max_merge_mm':1e-8/2**level}
            if mode!='stable':
                plan=plans[seed];plan_path=out/(name+'.plan')
                plan_path.write_text('QP_PION_PLAN_V1\n'+''.join(str(q)+' '+' '.join(format(x,'.17g') for x in [p['proper_ps'],*p['direction']])+'\n' for q,p in enumerate(plan)))
                c.update(plan=plan,plan_file=str(plan_path),plan_sha256=hashlib.sha256(plan_path.read_bytes()).hexdigest())
            cases.append(c)
        for level in (7,8):add(0,'stable',level)
        for seed in range(1,seeds+1):
            for level in ((7,8) if seed==1 else (8,)):
                for mode in ('isotropic','collinear_forward','collinear_backward'):add(seed,mode,level)
    base.update(suite='triplet_ensemble',scope='FORCED_PION_KINEMATIC_GRID_WITH_CONDITIONAL_DIRECTION_CONTROLS',
        cases=cases,points=POINTS,sampled_seeds=list(range(1,seeds+1)),
        nominal_stop_ps=max(c['nominal_stop_ps'] for c in cases),
        control_definition='same proper-time thresholds; pion parallel/antiparallel to instantaneous parent momentum in the production CM; exact two-body kinematics',
        population_definition='only isotropic level8 seeds are the pilot population; points are not PDF/cross-section weighted',
        conditional_control_not_isotropic=True,physical_survival_probability=None,
        adapter_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ASSETS.iterdir() if p.is_file()})
    base['driver_sha256']={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in ('af_triplet_recoil_ensemble.py','af_quirkpheno_weak_pilot.py','audit_quirkpheno_sampled_pion.py')}
    (out/'manifest.json').write_text(json.dumps(base,indent=2)+'\n')
    return base


def first_minimum(a,case):
    beta=case['pair_boost_z'];gamma=1/math.sqrt(1-beta*beta);pair=[]
    # Search a broad first-return region before exporting a narrow window.
    centre=case['predicted_first_return_ps']
    for q in (0,1):
        ix=(a['q']==q)&(a['t']*1000>=.65*centre)&(a['t']*1000<=1.35*centre)
        t=a['t'][ix]*1000;x=np.column_stack([a[k][ix] for k in ('x','y','z')])
        if len(t)<2:raise ValueError('missing first-return search history')
        cm_t=gamma*(t-beta*x[:,2]/.299792458)
        x[:,2]=gamma*(x[:,2]-beta*.299792458*t)
        keep=np.r_[np.diff(cm_t)>0,True];pair.append((cm_t[keep],x[keep]))
    minimum=chord_minimum(pair)
    lo=max(t[0] for t,x in pair);hi=min(t[-1] for t,x in pair)
    if not lo+.01*(hi-lo)<minimum['time_ps']<hi-.01*(hi-lo):raise ValueError('minimum lies at search-window edge')
    # Convert the two endpoint coordinates at common CM time back to lab.
    labs=[]
    for t,x in pair:
        z=np.interp(minimum['time_ps'],t,x[:,2])
        labs.append(gamma*(minimum['time_ps']+beta*z/.299792458))
    return minimum,labs


def package(out):
    import uproot
    from audit_quirkpheno_sampled_pion import audit
    manifest=json.loads((out/'manifest.json').read_text());checked=audit(out)
    (out/'local_audit.json').write_text(json.dumps(checked,indent=2)+'\n')
    if not checked['local_checks_pass']:raise RuntimeError('full local audit failed; retain ROOT files and return local_audit.json plus logs')
    windows={};groups={};control_groups={}
    for case in manifest['cases']:
        with uproot.open(out/(case['name']+'.root')) as root:a=root['steps'].arrays(library='np')
        minimum,labs=first_minimum(a,case)
        beta=case['pair_boost_z'];gamma=1/math.sqrt(1-beta*beta)
        decays=checked['cases'][case['name']]['decays']
        before=[d for d in decays if gamma*(d['time_ps']-beta*d['position_mm'][2]/.299792458)<=minimum['time_ps']]
        case['decays_by_first_return']=len(before);case['return_minimum_cm_chord']=minimum
        case['window_fit_status']='BOTH_NEUTRAL_AT_FIRST_RETURN' if len(before)==2 else 'NOT_BOTH_NEUTRAL_AT_FIRST_RETURN'
        keep=np.zeros(len(a['q']),bool)
        for q in (0,1):
            candidate=(a['q']==q)&(abs(a['t']*1000-labs[q])<=.5)
            if len(before)==2:
                td=next(d['time_ps'] for d in before if (d['parent_pdg']>0)==(q==0))
                candidate&=(a['t']*1000>td)&(abs(a['pdg'])==9900017)
            keep|=candidate
        if len(before)==2:
            spans=[]
            for q in (0,1):
                ix=keep&(a['q']==q);tt=a['t'][ix]*1000;zz=a['z'][ix]
                cm=gamma*(tt-beta*zz/.299792458)
                spans.append((float(cm.min()),float(cm.max())) if len(cm)>=2 else (math.inf,-math.inf))
            if not max(s[0] for s in spans)<minimum['time_ps']<min(s[1] for s in spans):
                case['window_fit_status']='BOTH_NEUTRAL_BUT_WINDOW_DOES_NOT_BRACKET_RETURN'
                keep=np.zeros(len(a['q']),bool)
                for q in (0,1):keep|=(a['q']==q)&(abs(a['t']*1000-labs[q])<=.5)
        buf=io.BytesIO();np.savez_compressed(buf,**{k:v[keep] for k,v in a.items()});windows[case['name']+'.npz']=buf.getvalue()
        if case['mode']=='isotropic' and case['level']==8:
            bins=groups.setdefault(case['point'],{'population_seeds':0,'decays_by_return':{'0':0,'1':0,'2':0}})
            bins['population_seeds']+=1;bins['decays_by_return'][str(len(before))]+=1
        if case['mode'].startswith('collinear') and case['level']==8:
            key=case['point']+'/'+case['mode']
            bins=control_groups.setdefault(key,{'conditional_control_seeds':0,'decays_by_return':{'0':0,'1':0,'2':0}})
            bins['conditional_control_seeds']+=1;bins['decays_by_return'][str(len(before))]+=1
    manifest['ensemble_summary']=groups
    manifest['conditional_control_summary']=control_groups
    manifest['compact_window_definition']='per endpoint +/-0.5 ps around a common production-CM closest-approach estimate; neutral-qualified windows clipped after decay'
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    target=out.parent/(out.name+'_audit.tar.gz')
    with tarfile.open(target,'w:gz') as tar:
        for name in ('manifest.json','local_audit.json','build.log','patch_check.log','patch_apply.log','checkout.log'):
            if (out/name).exists():tar.add(out/name,arcname=name)
        for case in manifest['cases']:
            for suffix in ('.mac','.card','.log','.plan'):
                p=out/(case['name']+suffix)
                if p.exists():tar.add(p,arcname=p.name)
        for name,data in windows.items():
            member=tarfile.TarInfo(name);member.size=len(data);tar.addfile(member,io.BytesIO(data))
    return target


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--seeds',type=int,default=8);p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--pack-only',action='store_true');p.add_argument('--case-timeout-seconds',type=int,default=120)
    p.add_argument('--workers',type=int,default=1,help='Independent one-CPU cases in parallel; use allocated node CPUs')
    p.add_argument('--duration-ps',type=float,help='Optional longer transport horizon')
    a=p.parse_args();out=a.output_dir.expanduser().resolve()
    if a.case_timeout_seconds<1 or a.workers<1:raise SystemExit('case timeout and workers must be positive')
    if a.duration_ps is not None and (not math.isfinite(a.duration_ps) or a.duration_ps<=0):raise SystemExit('positive finite duration required')
    if a.prepare_only and a.pack_only:raise SystemExit('choose prepare-only or pack-only')
    if a.pack_only:
        saved=json.loads((out/'manifest.json').read_text())
        if saved.get('postneutral_packaging',saved.get('nominal_stop_ps',0)>=1000):
            from package_triplet_postneutral import package as package_late
            print('Upload:',package_late(out))
        else:print('Upload:',package(out))
        return
    if not a.prepare_only:
        try:import uproot  # noqa: F401
        except ImportError:raise SystemExit('uproot is required for automatic local audit and compact packaging')
    state=prepare(out,a.seeds)
    state['parallel_case_workers']=a.workers
    state['postneutral_packaging']=a.duration_ps is not None
    if a.duration_ps is not None:
        for case in state['cases']:
            path=Path(case['macro']);macro=re.sub(r'^/qsim/maxTime .*$',f'/qsim/maxTime {a.duration_ps/1000:.17g} ns',path.read_text(),flags=re.M)
            path.write_text(macro);case['macro_sha256']=hashlib.sha256(macro.encode()).hexdigest();case['nominal_stop_ps']=a.duration_ps
        state['nominal_stop_ps']=a.duration_ps
    (out/'manifest.json').write_text(json.dumps(state,indent=2)+'\n')
    print('Manifest:',out/'manifest.json','Cases:',len(state['cases']),flush=True)
    if a.prepare_only:return
    execute(state,out,a.case_timeout_seconds,False,a.workers)
    if a.duration_ps is not None:
        from package_triplet_postneutral import package as package_late
        print('Upload:',package_late(out),flush=True)
    else:print('Upload:',package(out),flush=True)


if __name__=='__main__':main()
