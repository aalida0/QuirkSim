"""16 direct parallel vacuum runs measuring endpoint+string orbital charges."""
import argparse,hashlib,json,re,sys
from pathlib import Path
from af_validate_inclusive_weak import prepare as inclusive_prepare
from af_quirkpheno_weak_pilot import execute
from audit_native_noether import audit_and_pack


def prepare(out):
    state=inclusive_prepare(out)
    selected=[]
    for c in state['cases']:
        if c['mode'] not in ('stable','isotropic'):
            for suffix in ('.mac','.card','.weakplan'):
                p=out/(c['name']+suffix)
                if p.exists():p.unlink()
            continue
        macro=out/(c['name']+'.mac');text=macro.read_text()
        # Seed/debt is a finite initialization regulator. Refine it as well.
        seed=1e-10/(2**(c['level']-7))
        for key,value in {'/quirk/firstStringLength':f'{seed:.17g} mm','/qsim/maxTime':'.15 ns','/qsim/out/stepPrescale':'0'}.items():
            text,n=re.subn('^'+re.escape(key)+r' .*$',key+' '+value,text,flags=re.M)
            if n!=1:raise ValueError('missing or ambiguous macro command '+key)
        macro.write_text(text);c.update(macro_sha256=hashlib.sha256(macro.read_bytes()).hexdigest(),nominal_stop_ps=150.,first_string_vector_length_mm=seed);selected.append(c)
    state.update(cases=selected,noether_diagnostic=True,suite='native_noether',scope='SAMPLED_CAUSAL_CONTOUR_CLASSICAL_CONSERVATION',nominal_stop_ps=150.,first_string_vector_length_mm='per_case',numerical_validation='NOT_RUN')
    root=Path(__file__).resolve().parents[1]
    for path in ('scripts/af_validate_native_noether.py','scripts/audit_native_noether.py'):
        state['driver_sha256'][path]=hashlib.sha256((root/path).read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps(state,indent=2)+'\n');return state


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',required=True,type=Path);p.add_argument('--workers',type=int,default=4);p.add_argument('--case-timeout-seconds',type=int,default=600);p.add_argument('--prepare-only',action='store_true');p.add_argument('--audit-only',action='store_true');a=p.parse_args()
    if a.workers<1 or a.case_timeout_seconds<1:p.error('positive workers and timeout required')
    out=a.output_dir.expanduser().resolve()
    if a.audit_only:print(audit_and_pack(out));sys.exit(0)
    state=prepare(out);print('Manifest:',out/'manifest.json','Cases:',len(state['cases']),flush=True)
    if not a.prepare_only:
        try:execute(state,out,a.case_timeout_seconds,workers=a.workers)
        finally:audit_and_pack(out)
