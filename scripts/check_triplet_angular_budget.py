"""Apply bounded-step reachability to saved neutral endpoint fits."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.angular_radiation import reachability


def calculate(source):
    data=json.loads(gzip.decompress(source.read_bytes()))
    model=data['runs'][0]['path']['model'];lam=model['lambda_ev']*1e-9
    binding=model['captured_binding_energy_gev'];rows=[]
    for fit in data['fits']:
        if fit['level']!=8:continue
        # Even credit the final compact binding release before low L is reached:
        # an intentionally generous budget, not an accessible cooling path.
        energy=fit['fitted_energy_gev']+binding
        for floor in (.1,1.,10.):
            for step in (1.,10.,100.):
                rows.append({'case':fit['case'],'point':fit['point'],'mode':fit['mode'],
                    'seed':fit['seed'],'minimum_energy_in_lambda':floor,'max_angular_step':step,
                    'generous_available_energy_gev':energy,
                    **reachability(fit['fitted_L_over_hbar'],energy,floor*lam,step)})
    return {'scope':'CONDITIONAL_BOUNDED_STEP_ANGULAR_REACHABILITY',
        'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'archive_sha256':data['archive_sha256'],'rows':rows,
        'assumptions':['Hard positive emission-energy floor and hard angular-step ceiling are declared.',
            'All available energy plus final compact binding release is credited; no time or contact gate imposed.',
            'Unbiased bound additionally requires zero conditional mean of vector angular increments.',
            'Endpoint central L is assumed to be the evolving action; string angular momentum is omitted.',
            'Massless emission, coherent torque or large-angular-momentum radiation can violate these assumptions.',
            'Zero here means inaccessible low L in this class; it is not a physical annihilation probability.']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=calculate(a.source);a.output.write_text(json.dumps(r,indent=2)+'\n')
    print('reachability cards',len(r['rows']))
