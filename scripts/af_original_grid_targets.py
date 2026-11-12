"""Bounded original-grid AF handoff: three pion-open points; no Condor.

This adapter is explicitly pion-pair conditioned. Its proper clocks use the
full LO weak width, while its omitted-channel population is kept in the card.
No results are fabricated by --prepare-only.
"""
import argparse,hashlib,json,math,re,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.physics import benchmark
from quirksim.decays import direction
from af_quirkpheno_weak_pilot import prepare_fixed,execute,ASSETS

TARGETS=[(100.,2154.),(225.,4642.),(400.,2154.)]

def prepare(out,seeds=4):
    if not 1<=seeds<=8:raise ValueError('bounded pilot requires 1..8 seeds')
    base=prepare_fixed(out);templates={level:(out/f'stable_parent_level{level}.mac').read_text() for level in (7,8)}
    for case in base['cases']:
        for suffix in ('.mac','.card'):
            path=out/(case['name']+suffix)
            if path.exists():path.unlink()
    cases=[];point_cards=[]
    for mass,lam in TARGETS:
        b=benchmark('triplet',mass);parent=b['charged_mass_gev'];br=b['widths_gev']['pion']/sum(b['widths_gev'].values())
        if br<=0:raise ValueError('C++ pion adapter cannot represent a closed pion channel')
        energy=parent/math.sqrt(.84);px=.4*energy;boost=.25;gamma=1/math.sqrt(1-boost*boost);pz=gamma*boost*energy
        point=f'm{mass:g}_lambda{lam:g}';point_cards.append({'name':point,'mass_gev':mass,'lambda_ev':lam,'weak':b,'pion_pair_fraction':br*br,'omitted_channel_pair_fraction':1-br*br})
        plans={}
        for seed in range(1,seeds+1):
            rng=np.random.default_rng(seed);thresholds=rng.exponential(b['lifetime_ps'],2)
            plans[seed]=[{'proper_ps':float(thresholds[i]),'direction':direction(rng).tolist()} for i in (0,1)]
        for field in (0.,2.):
            for mode in ('stable','isotropic','collinear_forward','collinear_backward'):
                for seed in ([0] if mode=='stable' else range(1,seeds+1)):
                    for level in ((7,8) if seed in (0,1) else (8,)):
                        name=f'{point}_B{field:g}_{mode}_seed{seed}_level{level}'
                        card=out/(name+'.card');card.write_text(f'world 100 100 G4_Galactic\nfield {field:g} 0 100 0 100\n')
                        macro=templates[level]
                        changes={'/quirk/mass':f'{parent:.17g} GeV','/quirk/lambda':f'{lam:.17g} eV',
                            '/qsim/gun/p1':f'{px:.17g} 0 {pz:.17g} GeV','/qsim/gun/p2':f'{-px:.17g} 0 {pz:.17g} GeV',
                            '/qsim/maxTime':'4 ns','/qsim/card':str(card),'/qsim/out/file':str(out/(name+'.root'))}
                        for key,value in changes.items():macro=re.sub('^'+re.escape(key)+r' .*$',key+' '+value,macro,flags=re.M)
                        path=out/(name+'.mac');path.write_text(macro)
                        case={'name':name,'point':point,'mode':mode,'level':level,'seed':seed,'macro':str(path),
                            'macro_sha256':hashlib.sha256(macro.encode()).hexdigest(),'weak_enabled':mode!='stable','endpoint':-1,
                            'clock_mode':'localized','neutral_mass_gev':mass,'mass_gev':parent,'lambda_ev':lam,
                            'field_tesla':field,'pair_boost_z':boost,'beta_each_cm':.4,'nominal_stop_ps':4000.,
                            'direction_mode':'collinear_cm' if mode.startswith('collinear') else 'planned',
                            'collinear_sign':-1 if mode=='collinear_backward' else 1,
                            'predicted_first_return_ps':2*px*6.582119569e-13/(lam*1e-9)**2*gamma}
                        if mode!='stable':
                            case['plan']=plans[seed];plan_path=out/(name+'.plan')
                            plan_path.write_text('QP_PION_PLAN_V1\n'+''.join(str(i)+' '+' '.join(f'{v:.17g}' for v in [p['proper_ps'],*p['direction']])+'\n' for i,p in enumerate(plans[seed])))
                            case.update(plan_file=str(plan_path),plan_sha256=hashlib.sha256(plan_path.read_bytes()).hexdigest())
                        cases.append(case)
    base.update(suite='triplet_ensemble',scope='ORIGINAL_GRID_TARGETED_PION_PAIR_CONDITIONED_TRANSPORT_NOT_ANNIHILATION',
        cases=cases,points=point_cards,nominal_stop_ps=4000.,pair_boost_z=.25,
        mass_gev=225.1571702763212,neutral_mass_gev=225.,
        mass_note='Legacy globals are fallbacks only; every case carries its actual parent/neutral mass and tension.',
        full_weak_channels_supported=False,physical_survival_probability=None,numerical_validation='NOT_RUN',
        driver_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ('af_original_grid_targets.py','af_quirkpheno_weak_pilot.py','audit_quirkpheno_sampled_pion.py','package_triplet_postneutral.py')},
        adapter_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ASSETS.iterdir() if p.is_file()})
    (out/'manifest.json').write_text(json.dumps(base,indent=2)+'\n');return base

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--seeds',type=int,default=4);p.add_argument('--workers',type=int,default=4)
    p.add_argument('--case-timeout-seconds',type=int,default=600);p.add_argument('--prepare-only',action='store_true');a=p.parse_args()
    out=a.output_dir.expanduser().resolve();state=prepare(out,a.seeds);print('Manifest:',out/'manifest.json','Cases:',len(state['cases']),flush=True)
    if not a.prepare_only:
        execute(state,out,a.case_timeout_seconds,workers=a.workers)
        from package_triplet_postneutral import package
        print('Compact archive:',package(out))
