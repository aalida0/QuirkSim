"""Execute pinned author formulas as a scoped external normalization check.

Pass the unmodified sommerfeld.py from the pinned public repository. The
legacy Python-2 CLI is excluded; its function definitions are not rewritten.
This is an SU(3), color/spin-averaged l=0 control, not the IC-singlet model.
"""
import argparse
import hashlib
import json
from pathlib import Path
import mpmath as mp

PIN='f0c97f32d10bb0f23f791765e652426117211955'
BLOB='982db5d710283acced8e4d5b580784fd6807d93b'


def check(source):
    raw=source.read_bytes()
    observed=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    if observed!=BLOB:raise ValueError('author source does not match pinned git blob')
    namespace={}
    exec(compile(raw.decode().split('# argument parser')[0],str(source),'exec'),namespace)
    cases=[]
    with mp.workdps(40):
        mass=mp.mpf(225);alpha=mp.mpf('.04')
        # The original CLI sets this global; generated functions use it despite
        # their similarly named formal argument. Reproduce the original context.
        namespace['alphasommerfeld']=alpha
        for beta in map(mp.mpf,('0.00001','.05','.4')):
            plain=namespace['get_xsec']('fftogg',3,0,False,mass,beta,alpha,alpha)
            corrected=namespace['get_xsec']('fftogg',3,0,True,mass,beta,alpha,alpha)
            reference=7*mp.pi*alpha**2/(27*mass**2)
            # Coulomb relative momentum p/m=gamma*beta in the author formulas.
            def factor(strength):
                z=mp.pi*strength*mp.sqrt(1-beta*beta)/beta
                return z/(-mp.expm1(-z))
            expected=mp.mpf(2)/7*factor(mp.mpf(4)/3*alpha)+mp.mpf(5)/7*factor(-alpha/6)
            bare_error=abs(plain*2*beta/reference-1)
            corrected_error=abs(corrected/plain/expected-1)
            if max(bare_error,corrected_error)>mp.mpf('1e-12'):
                raise RuntimeError('published normalization control failed')
            cases.append({'beta':float(beta),'sigma_v_moller_gev_minus2':float(plain*2*beta),
                          'sommerfeld_ratio':float(corrected/plain),
                          'bare_relative_error':float(bare_error),'coulomb_relative_error':float(corrected_error)})
    return {'status':'PASS_SCOPED_AUTHOR_CODE_CONTROL','repository':'https://github.com/MDT-Maikel/Sommerfeld',
            'commit':PIN,'git_blob':observed,'source_sha256':hashlib.sha256(raw).hexdigest(),
            'license':'CC-BY-4.0 per upstream README',
            'coverage':'SU3_FUNDAMENTAL_COLOR_AND_SPIN_AVERAGED_L0_ONLY',
            'source_processing':'execute definitions before legacy CLI, reproduce CLI global alphasommerfeld',
            'not_validated':['SU2 general-N notebook execution','IC-singlet partial-wave projection',
                             'capture, confinement or nonperturbative emission'],
            'cases':cases}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=check(a.source);a.output.write_text(json.dumps(r,indent=2)+'\n');print(r['status'])
