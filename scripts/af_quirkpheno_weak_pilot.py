"""Bounded deterministic one-pion decay integration pilot; no survival claim."""
import argparse
import hashlib
import json
import math
import random
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor,as_completed
from af_quirkpheno_transport import prepare as stable_prepare, run_logged, LCG, PIN, URL

M_NEUTRAL = 225.
M_PARENT = 225.1571702763212
ASSETS = Path(__file__).resolve().parents[1] / 'research' / 'quirkpheno_recoil'

def prepare_fixed(out, clock_mode="localized"):
    state = stable_prepare(out, 'precision')
    templates = {c['level']: (out / (c['name']+'.mac')).read_text()
                 for c in state['cases'] if c['field_tesla']==2 and c['level'] in (7,8)}
    for c in state['cases']:
        (out / (c['name']+'.mac')).unlink()
        (out / (c['name']+'.card')).unlink()
    cases=[]
    energy=M_PARENT/math.sqrt(1-.4**2)
    px=.4*energy; pz=.25*energy/math.sqrt(1-.25**2)
    for mode,endpoint,sign in [('stable_parent',0,1),('pion_yplus_end0',0,1),
                                ('pion_yminus_end0',0,-1),('pion_yplus_end1',1,1)]:
        for level in (7,8):
            name=f'{mode}_level{level}'
            card=out/(name+'.card'); card.write_text('world 100 100 G4_Galactic\nfield 2 0 100 0 100\n')
            macro=templates[level]
            macro=re.sub(r'^/quirk/mass .*$',f'/quirk/mass {M_PARENT:.17g} GeV',macro,flags=re.M)
            macro=re.sub(r'^/qsim/gun/p1 .*$',f'/qsim/gun/p1 {px:.17g} 0 {pz:.17g} GeV',macro,flags=re.M)
            macro=re.sub(r'^/qsim/gun/p2 .*$',f'/qsim/gun/p2 {-px:.17g} 0 {pz:.17g} GeV',macro,flags=re.M)
            macro=re.sub(r'^/qsim/card .*$',f'/qsim/card {card}',macro,flags=re.M)
            macro=re.sub(r'^/qsim/out/file .*$',f'/qsim/out/file {out/(name+".root")}',macro,flags=re.M)
            path=out/(name+'.mac');path.write_text(macro)
            cases.append({'name':name,'mode':mode,'level':level,'macro':str(path),
                          'macro_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                          'field_tesla':2,'weak_enabled':mode!='stable_parent',
                          'endpoint':endpoint,'pion_rest_y_sign':sign,'proper_threshold_ps':50.,
                          'max_boost':.01/2**level,'max_merge_mm':1e-8/2**level,'clock_mode':clock_mode})
    if clock_mode=='localized':
        for original_mode in ('stable_parent','pion_yplus_end0'):
            template=next(c for c in cases if c['mode']==original_mode and c['level']==8)
            name='clock_legacy_'+template['name']
            card=out/(name+'.card');card.write_text((out/(template['name']+'.card')).read_text())
            macro=(out/(template['name']+'.mac')).read_text().replace(template['name'],name)
            path=out/(name+'.mac');path.write_text(macro)
            case=dict(template,name=name,mode='clock_legacy_'+original_mode,clock_mode='legacy',macro=str(path),
                      macro_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            cases.append(case)
    state.update(suite='deterministic_weak_handoff',scope='EXPERIMENTAL_SINGLE_PION_HANDOFF_NOT_LIFETIME_ENSEMBLE',
                 mass_gev=M_PARENT,neutral_mass_gev=M_NEUTRAL,delta_mev=(M_PARENT-M_NEUTRAL)*1000,
                 mass_convention='existing QuirkSim LO triplet benchmark; parent mass is loop split',
                 decay_clock='localized Gauss3 m/E clock and event root; two same-input legacy controls' if clock_mode=='localized' else 'legacy first post-step threshold',
                 default_clock_mode=clock_mode,proper_root_tolerance_ps=1e-9 if clock_mode=='localized' else None,
                 modes_are_physical_angles=True,causal_recoil_validation='PENDING',cases=cases)
    state['adapter_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in sorted(ASSETS.iterdir()) if p.is_file()}
    return state

def prepare(out, clock_mode="localized", sampled_pion_seeds=0):
    if not sampled_pion_seeds:
        return prepare_fixed(out,clock_mode)
    if clock_mode!='localized' or not 1<=sampled_pion_seeds<=32:
        raise ValueError('Sampled pion pilot requires localized clock and 1-32 complete seeds')
    state=prepare_fixed(out,clock_mode)
    templates={c['level']:(out/(c['name']+'.mac')).read_text() for c in state['cases'] if c['mode']=='stable_parent'}
    for c in state['cases']:
        (out/(c['name']+'.mac')).unlink();(out/(c['name']+'.card')).unlink()
    rate_file=ASSETS/'triplet_pion_rate_reference.json'
    rates=json.loads(rate_file.read_text())
    if abs(rates['charged_mass_gev']-M_PARENT)>1e-12 or rates['neutral_mass_gev']!=M_NEUTRAL:
        raise ValueError('Frozen rate reference and transport masses disagree')
    cases=[]
    def add(name,mode,level,plan=None,seed=None):
        card=out/(name+'.card');card.write_text('world 100 100 G4_Galactic\nfield 2 0 100 0 100\n')
        macro=templates[level].replace('stable_parent_level'+str(level),name)
        path=out/(name+'.mac');path.write_text(macro)
        case={'name':name,'mode':mode,'level':level,'macro':str(path),
              'macro_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'field_tesla':2,
              'weak_enabled':plan is not None,'endpoint':-1 if plan else 0,'seed':seed,
              'clock_mode':'localized','max_boost':.01/2**level,'max_merge_mm':1e-8/2**level}
        if plan:
            plan_path=out/(name+'.plan')
            text='QP_PION_PLAN_V1\n'+''.join(str(i)+' '+' '.join(format(x,'.17g') for x in [p['proper_ps'],*p['direction']])+'\n' for i,p in enumerate(plan))
            plan_path.write_text(text)
            case.update(plan_file=str(plan_path),plan_sha256=hashlib.sha256(plan_path.read_bytes()).hexdigest(),plan=plan)
        cases.append(case)
    for level in (7,8):add('sampled_stable_level'+str(level),'sampled_stable',level)
    for seed in range(1,sampled_pion_seeds+1):
        rng=random.Random(seed);plan=[]
        for endpoint in (0,1):
            proper=rng.expovariate(1/rates['lifetime_ps']);costheta=rng.uniform(-1,1);phi=rng.uniform(0,2*math.pi)
            st=math.sqrt(max(0.,1-costheta*costheta))
            plan.append({'proper_ps':proper,'direction':[st*math.cos(phi),st*math.sin(phi),costheta]})
        for level in ((7,8) if seed==1 else (8,)):
            for opposite in (False,True):
                mode=('antipodal' if opposite else 'sampled')+'_seed'+str(seed)
                selected=[{'proper_ps':p['proper_ps'],'direction':[(-1 if opposite else 1)*x for x in p['direction']]} for p in plan]
                add(mode+'_level'+str(level),mode,level,selected,seed)
    state.update(suite='sampled_pion',scope='FORCED_PION_CHANNEL_FIXED_KINEMATICS_FINITE_WINDOW_NO_SURVIVAL',
                 sampled_seeds=list(range(1,sampled_pion_seeds+1)),cases=cases,
                 pion_only_conditional=True,inclusive_channel_sampling=False,
                 lifetime_mean_ps=rates['lifetime_ps'],reference_pion_branching=rates['pion_branching'],
                 rate_reference_sha256=hashlib.sha256(rate_file.read_bytes()).hexdigest(),
                 decay_clock='independent exponential proper thresholds for both endpoints; localized roots',
                 random_generator='Python random.Random(seed); per endpoint expovariate, cos(theta), phi')
    return state

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',required=True)
    p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--clock-mode',choices=('localized','legacy'),default='localized')
    p.add_argument('--sampled-pion-seeds',type=int,default=0,help='Optional complete forced-pion seed block; 8 seeds gives 20 cases')
    p.add_argument('--duration-ps',type=float,default=130.)
    p.add_argument('--case-timeout-seconds',type=int,default=120)
    p.add_argument('--audit-close-return',action='store_true',help='Audit full step records after each run; requires uproot in this Python environment')
    a=p.parse_args();out=Path(a.output_dir).expanduser().resolve();state=prepare(out,a.clock_mode,a.sampled_pion_seeds)
    if not math.isfinite(a.duration_ps) or a.duration_ps<=5 or a.case_timeout_seconds<1:
        raise SystemExit('duration must be finite and >5 ps and timeout must be positive')
    if a.audit_close_return:
        try: import uproot  # noqa: F401
        except ImportError: raise SystemExit('--audit-close-return requires uproot in this Python environment')
    for case in state['cases']:
        macro=Path(case['macro'])
        content=re.sub(r'^/qsim/maxTime .*$',f'/qsim/maxTime {a.duration_ps/1000:.17g} ns',macro.read_text(),flags=re.M)
        macro.write_text(content);case['macro_sha256']=hashlib.sha256(content.encode()).hexdigest()
    state['nominal_stop_ps']=a.duration_ps;state['timeout_seconds_per_case']=a.case_timeout_seconds
    state['close_return_audit_enabled']=a.audit_close_return
    def save(): (out/'manifest.json').write_text(json.dumps(state,indent=2)+'\n')
    save();print(f'Manifest: {out/"manifest.json"}',flush=True)
    if a.prepare_only:return
    execute(state,out,a.case_timeout_seconds,a.audit_close_return)

def execute(state,out,case_timeout_seconds=120,audit_close_return=False,workers=1):
    """Build the pinned adapter once and run an already prepared manifest."""
    if type(workers) is not int or workers<1:raise ValueError('workers must be a positive integer')
    state['parallel_case_workers']=workers
    def save(): (out/'manifest.json').write_text(json.dumps(state,indent=2)+'\n')
    if not LCG.is_file():raise SystemExit(f'Missing Geant4 setup: {LCG}')
    source=out/'source'
    def checked(cmd,cwd,log,seconds):
        r=run_logged(cmd,cwd,out/log,seconds)
        state[log]=r;save()
        if r['returncode']!=0:
            print('\n'.join((out/log).read_text(errors='replace').splitlines()[-35:]),file=sys.stderr)
            raise SystemExit(f'Failed: {log}; retain output directory')
    checked(['git','clone',URL,str(source)],out,'clone.log',120)
    checked(['git','checkout','--detach',PIN],source,'checkout.log',60)
    state['observed_source_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip()
    if state['observed_source_commit']!=PIN:
        raise SystemExit('Unexpected upstream revision; refusing patch')
    preimages=json.loads((ASSETS/'upstream_preimages.json').read_text())
    state['upstream_preimage_checks']={}
    for relative,expected in preimages['files'].items():
        actual=hashlib.sha256((source/relative).read_bytes()).hexdigest()
        state['upstream_preimage_checks'][relative]=actual==expected['sha256']
    save()
    if not all(state['upstream_preimage_checks'].values()):
        raise SystemExit('Upstream source bytes differ from pinned preimages; refusing patch')
    patch=ASSETS/'weak_handoff.patch' 
    checked(['git','apply','--check',str(patch)],source,'patch_check.log',60)
    checked(['git','apply',str(patch)],source,'patch_apply.log',60)
    for name in ('EndpointHandoff.h','TwoBodyKinematics.h','ProperTimeClock.h','SampledPionPlan.h','DirectionControl.h','GeneralWeakPlan.h','QuirkWeakProcess.h'):
        shutil.copyfile(ASSETS/name,source/'quirkcore/include/quirkcore'/name)
    shutil.copyfile(ASSETS/'QuirkWeakProcess.cc',source/'quirkcore/src/QuirkWeakProcess.cc')
    build=out/'build'
    def lcg(cmd, timeout_seconds=600):
        return ['timeout','--kill-after=10s',str(timeout_seconds)+'s','bash','-c',
                'source '+shlex.quote(str(LCG))+' && '+cmd]
    cmd='cmake -S '+shlex.quote(str(source))+' -B '+shlex.quote(str(build))
    cmd+=' -DQSIM_WITH_HEPMC3=OFF -DQSIM_WITH_PYTHIA8=OFF && cmake --build '+shlex.quote(str(build))+' -j1'
    checked(lcg(cmd),source,'build.log',620)
    state['build']=state['build.log'];save()
    exe=build/'quirksim/quirksim'
    def run_case(case):
        env={'QP_WEAK_ENABLED':int(case['weak_enabled']),'QP_WEAK_ENDPOINT':case['endpoint'],
             'QP_WEAK_TAU_PS':50,'QP_WEAK_CLOCK_MODE':case['clock_mode'],'QP_NEUTRAL_MASS_GEV':case.get('neutral_mass_gev',M_NEUTRAL),'QP_WEAK_DIR_SIGN':case.get('pion_rest_y_sign',1),'QP_WEAK_PLAN':case.get('plan_file','')}
        env['QP_WEAK_GENERAL_PLAN']=case.get('general_plan_file','')
        env.update(QP_WEAK_DIRECTION_MODE=case.get('direction_mode','planned'),QP_WEAK_COLLINEAR_SIGN=case.get('collinear_sign',1))
        env['QP_WEAK_CONTROL_BOOST_Z']=case.get('pair_boost_z',0.)
        command='env '+' '.join(shlex.quote(k+'='+str(v)) for k,v in env.items())
        command+=' timeout --kill-after=5s '+str(case_timeout_seconds)+'s '+shlex.quote(str(exe))
        command+=' --mt -t 1 -n 1 -s 1 -m '+shlex.quote(case['macro'])
        execution=run_logged(lcg(command,case_timeout_seconds+15),source,out/(case['name']+'.log'),case_timeout_seconds+25)
        return case,execution,(out/(case['name']+'.root')).is_file()
    pool=ThreadPoolExecutor(max_workers=workers)
    pending=[pool.submit(run_case,case) for case in state['cases']]
    try:
      for future in as_completed(pending):
        case,case['execution'],case['root_exists']=future.result()
        print('Finished '+case['name'],flush=True)
        save()
        if case['execution']['returncode']!=0 or not case['root_exists']:
            state['numerical_validation']='STOPPED_AFTER_FAILED_CASE';save()
            raise SystemExit('Failed/incomplete: '+case['name']+'; retain and return logs')
        if audit_close_return:
            from audit_quirkpheno_close_return import audit
            try: case['close_return']=audit(out/(case['name']+'.root'))
            except (ValueError, KeyError, OSError) as exc:
                case['close_return']={'status':'INVALID_INPUT','error':str(exc)}
                state['numerical_validation']='STOPPED_AFTER_INVALID_CLOSE_RETURN_INPUT';save()
                raise SystemExit('Close-return audit failed: '+case['name']) from exc
            save()
    finally:
        pool.shutdown(wait=True,cancel_futures=True)
    state['numerical_validation']='PENDING_WEAK_ROOT_AND_HANDOFF_AUDIT';save()
    failed=[c['name'] for c in state['cases'] if c['execution']['returncode']!=0 or not c['root_exists']]
    if failed:raise SystemExit('Failed/incomplete: '+', '.join(failed)+'; return logs')
    print('Transport complete; preparing compact local audit.' if state.get('suite')=='triplet_ensemble' else
          'Finished. Upload manifest, macros, logs and ROOT files for handoff audit.')
if __name__=='__main__':main()
