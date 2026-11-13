"""40-case paired native weak-handoff validation, direct parallel AF workers.

Forced channels and proper clocks test the adapter; this is NOT a lifetime
ensemble, production sample, or annihilation calculation. ROOT stays on AF.
"""
import argparse,hashlib,json,math,re,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.physics import benchmark
from quirksim.weak_plan import make_plan,write_plan
from af_quirkpheno_weak_pilot import prepare_fixed,execute,ASSETS


def prepare(out):
    state=prepare_fixed(out);templates={k:(out/f'stable_parent_level{k}.mac').read_text() for k in (7,8)}
    for c in state['cases']:
        for suffix in ('.mac','.card'):
            p=out/(c['name']+suffix)
            if p.exists():p.unlink()
    b=benchmark('triplet',225.);M=b['charged_mass_gev'];en=M/math.sqrt(.84);boost=.25;gam=1/math.sqrt(1-boost*boost);pz=gam*boost*en
    cases=[]
    for channel in ('stable','pion','electron_neutrino','muon_neutrino'):
      plan=None if channel=='stable' else make_plan(225.,31,channel=channel,proper_ps=(50.,100.))
      for field in (0.,2.):
       for mode in (('stable',) if channel=='stable' else ('isotropic','collinear_forward','collinear_backward')):
        for level in (7,8):
         name=f'{channel}_B{field:g}_{mode}_level{level}';card=out/(name+'.card');card.write_text(f'world 100 100 G4_Galactic\nfield {field:g} 0 100 0 100\n')
         text=templates[level]
         changes={'/quirk/mass':f'{M:.17g} GeV','/quirk/lambda':'4642 eV','/qsim/gun/p1':f'{.4*en:.17g} 0 {pz:.17g} GeV','/qsim/gun/p2':f'{-.4*en:.17g} 0 {pz:.17g} GeV','/qsim/maxTime':'.5 ns','/qsim/card':str(card),'/qsim/out/file':str(out/(name+'.root'))}
         for k,v in changes.items():text=re.sub('^'+re.escape(k)+r' .*$',k+' '+v,text,flags=re.M)
         macro=out/(name+'.mac');macro.write_text(text)
         case=dict(name=name,point='m225_lambda4642',mode=mode,channel=channel,level=level,seed=31,weak_enabled=plan is not None,endpoint=-1,clock_mode='localized',neutral_mass_gev=225.,mass_gev=M,lambda_ev=4642.,field_tesla=field,pair_boost_z=boost,beta_each_cm=.4,nominal_stop_ps=500.,macro=str(macro),macro_sha256=hashlib.sha256(macro.read_bytes()).hexdigest(),direction_mode='collinear_cm' if mode.startswith('collinear') else 'planned',collinear_sign=-1 if mode=='collinear_backward' else 1)
         if plan:
            path=out/(name+'.weakplan');write_plan(path,plan);case.update(general_plan_file=str(path),general_plan_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),plan=plan)
         cases.append(case)
    state.update(scope='FORCED_CHANNEL_INCLUSIVE_VERTEX_AND_STRING_HANDOFF_VALIDATION_NOT_POPULATION',suite='inclusive_weak_validation',cases=cases,nominal_stop_ps=500.,lambda_ev=4642.,weak_reference=b,full_weak_channels_supported=True,numerical_validation='NOT_RUN',physical_survival_probability=None,adapter_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ASSETS.iterdir() if p.is_file()},validation_clocks_ps=[50.,100.],standard_product_mass_convention='Plan masses explicitly assigned to Geant4 dynamic particles; table agreement required within 1 keV')
    root=Path(__file__).resolve().parents[1]
    state['driver_sha256']={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ('scripts/af_validate_inclusive_weak.py','scripts/af_quirkpheno_weak_pilot.py','scripts/audit_inclusive_weak.py','quirksim/weak_plan.py','quirksim/decays.py','quirksim/physics.py')}
    (out/'manifest.json').write_text(json.dumps(state,indent=2)+'\n');return state


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',required=True,type=Path);p.add_argument('--workers',type=int,default=4);p.add_argument('--case-timeout-seconds',type=int,default=600);p.add_argument('--prepare-only',action='store_true');p.add_argument('--audit-only',action='store_true');a=p.parse_args()
    if a.workers<1 or a.case_timeout_seconds<1:p.error('positive worker count and timeout required')
    out=a.output_dir.resolve()
    if a.audit_only:
        from audit_inclusive_weak import audit_and_pack
        print('Compact audit:',audit_and_pack(out));sys.exit(0)
    if not a.prepare_only:
        try:import uproot
        except ImportError:raise SystemExit('Install uproot in the active Python environment before starting: python -m pip install uproot')
    state=prepare(out);print('Manifest:',out/'manifest.json','Cases:',len(state['cases']),flush=True)
    if not a.prepare_only:
        from audit_inclusive_weak import audit_and_pack
        try:execute(state,out,a.case_timeout_seconds,workers=a.workers)
        except SystemExit:
            print('Partial compact audit:',audit_and_pack(out),flush=True)
            raise
        print('Compact audit:',audit_and_pack(out),flush=True)
