"""Bounded deterministic one-pion decay integration pilot; no survival claim."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from af_quirkpheno_transport import prepare as stable_prepare, run_logged, LCG, PIN, URL

M_NEUTRAL = 225.
M_PARENT = 225.1571702763212
ASSETS = Path(__file__).resolve().parents[1] / 'research' / 'quirkpheno_recoil'

def prepare(out, clock_mode="localized"):
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

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',required=True)
    p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--clock-mode',choices=('localized','legacy'),default='localized')
    a=p.parse_args();out=Path(a.output_dir).expanduser().resolve();state=prepare(out,a.clock_mode)
    def save(): (out/'manifest.json').write_text(json.dumps(state,indent=2)+'\n')
    save();print(f'Manifest: {out/"manifest.json"}',flush=True)
    if a.prepare_only:return
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
    for name in ('EndpointHandoff.h','TwoBodyKinematics.h','ProperTimeClock.h','QuirkWeakProcess.h'):
        shutil.copyfile(ASSETS/name,source/'quirkcore/include/quirkcore'/name)
    shutil.copyfile(ASSETS/'QuirkWeakProcess.cc',source/'quirkcore/src/QuirkWeakProcess.cc')
    build=out/'build'
    def lcg(cmd):
        return ['timeout','--kill-after=10s','600s','bash','-c',
                'source '+shlex.quote(str(LCG))+' && '+cmd]
    cmd='cmake -S '+shlex.quote(str(source))+' -B '+shlex.quote(str(build))
    cmd+=' -DQSIM_WITH_HEPMC3=OFF -DQSIM_WITH_PYTHIA8=OFF && cmake --build '+shlex.quote(str(build))+' -j1'
    checked(lcg(cmd),source,'build.log',620)
    state['build']=state['build.log'];save()
    exe=build/'quirksim/quirksim'
    for case in state['cases']:
        print('Running '+case['name'],flush=True)
        env={'QP_WEAK_ENABLED':int(case['weak_enabled']),'QP_WEAK_ENDPOINT':case['endpoint'],
             'QP_WEAK_TAU_PS':50,'QP_WEAK_CLOCK_MODE':case['clock_mode'],'QP_NEUTRAL_MASS_GEV':M_NEUTRAL,'QP_WEAK_DIR_SIGN':case['pion_rest_y_sign']}
        command='env '+' '.join(shlex.quote(k+'='+str(v)) for k,v in env.items())
        command+=' timeout --kill-after=5s 120s '+shlex.quote(str(exe))
        command+=' --mt -t 1 -n 1 -s 1 -m '+shlex.quote(case['macro'])
        case['execution']=run_logged(lcg(command),source,out/(case['name']+'.log'),145)
        case['root_exists']=(out/(case['name']+'.root')).is_file()
        save()
        if case['execution']['returncode']!=0 or not case['root_exists']:
            state['numerical_validation']='STOPPED_AFTER_FAILED_CASE';save()
            raise SystemExit('Failed/incomplete: '+case['name']+'; retain and return logs')
    state['numerical_validation']='PENDING_WEAK_ROOT_AND_HANDOFF_AUDIT';save()
    failed=[c['name'] for c in state['cases'] if c['execution']['returncode']!=0 or not c['root_exists']]
    if failed:raise SystemExit('Failed/incomplete: '+', '.join(failed)+'; return logs')
    print('Finished. Upload manifest, macros, logs and ROOT files for handoff audit.')
if __name__=='__main__':main()
