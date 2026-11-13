"""Small inclusive LO decay-time/channel ensemble with paired physical controls.

Fixed production point, vacuum 0/2 T; no radiation/capture or pp-population claim.
"""
import argparse,copy,hashlib,json,math,re,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.weak_plan import make_plan,write_plan
from af_validate_native_noether import prepare as prepare_reference
from af_quirkpheno_weak_pilot import execute
from audit_physical_weak import audit_and_pack


def prepare(out,seeds=8,duration=1000.):
    if type(seeds) is not int or not 1<=seeds<=64 or not math.isfinite(duration) or duration<=150:raise ValueError('require 1..64 seeds and finite duration >150 ps')
    state=prepare_reference(out);templates={}
    cases=[]
    for c in state['cases']:
        if c['channel']=='stable':
            c['seed']=None;templates[(c['field_tesla'],c['level'])]=copy.deepcopy(c);cases.append(c)
        else:
            for suffix in ('.mac','.card','.weakplan'):
                p=out/(c['name']+suffix)
                if p.exists():p.unlink()
    for seed in range(1,seeds+1):
      plan=make_plan(225.,seed)  # No forced channel or proper clock.
      for field in (0.,2.):
       for mode in ('recoil','collinear_forward','collinear_backward'):
        for level in ((7,8) if seed==1 else (8,)):
            template=templates[(field,level)];name=f'physical_seed{seed}_B{field:g}_{mode}_level{level}'
            c=copy.deepcopy(template);c.update(name=name,seed=seed,mode=mode,channel='inclusive',weak_enabled=True,endpoint=-1,direction_mode='planned' if mode=='recoil' else 'collinear_cm',collinear_sign=-1 if mode=='collinear_backward' else 1,plan=plan)
            for suffix in ('.mac','.card'):
                p=out/(name+suffix);p.write_text((out/(template['name']+suffix)).read_text().replace(template['name'],name))
            c['macro']=str(out/(name+'.mac'));path=out/(name+'.weakplan');write_plan(path,plan)
            c.update(general_plan_file=str(path),general_plan_sha256=hashlib.sha256(path.read_bytes()).hexdigest());cases.append(c)
    for c in cases:
        p=Path(c['macro']);text,n=re.subn(r'^/qsim/maxTime .*$',f'/qsim/maxTime {duration/1000:.17g} ns',p.read_text(),flags=re.M)
        if n!=1:raise ValueError('missing stop command')
        p.write_text(text);c.update(macro_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),nominal_stop_ps=duration)
    state.update(cases=cases,physical_weak_ensemble=True,suite='physical_weak_ensemble',scope='INCLUSIVE_LO_WEAK_FIXED_PRODUCTION_FINITE_WINDOW_NO_RADIATION_NO_CAPTURE',nominal_stop_ps=duration,sampled_seeds=list(range(1,seeds+1)),validation_clocks_ps=None,decay_clock='independent exponential proper times and inclusive channels from benchmark LO rates',independent_events_per_field=seeds,paired_controls_are_not_independent_events=True,numerical_validation='NOT_RUN')
    state['conservation_preflight_cases']=[c['name'] for c in cases if c['seed']==1 and c['field_tesla']==0 and c['mode']=='recoil']
    root=Path(__file__).resolve().parents[1]
    for path in ('scripts/af_physical_weak_ensemble.py','scripts/audit_physical_weak.py'):
        state['driver_sha256'][path]=hashlib.sha256((root/path).read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps(state,indent=2)+'\n');return state


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',required=True,type=Path);p.add_argument('--workers',type=int,default=4);p.add_argument('--seeds',type=int,default=8);p.add_argument('--duration-ps',type=float,default=1000.);p.add_argument('--case-timeout-seconds',type=int,default=1200);p.add_argument('--prepare-only',action='store_true');p.add_argument('--audit-only',action='store_true');a=p.parse_args()
    if a.workers<1 or a.case_timeout_seconds<1:p.error('positive workers and timeout required')
    out=a.output_dir.expanduser().resolve()
    if a.audit_only:print(audit_and_pack(out));sys.exit(0)
    state=prepare(out,a.seeds,a.duration_ps);print('Manifest:',out/'manifest.json','Cases:',len(state['cases']),flush=True)
    if not a.prepare_only:
        try:execute(state,out,a.case_timeout_seconds,workers=a.workers)
        finally:audit_and_pack(out)
