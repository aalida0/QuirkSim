"""Audit small NPZ return-window archives without interpolating for exclusion."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.close_return import speed_cone_window


def chord_minimum(pair):
    lo=max(t[0] for t,x in pair);hi=min(t[-1] for t,x in pair)
    t=np.unique(np.concatenate([v[0] for v in pair]));t=t[(t>=lo)&(t<=hi)]
    x=[np.column_stack([np.interp(t,v[0],v[1][:,j]) for j in range(3)]) for v in pair]
    r=x[0]-x[1];dr=np.diff(r,axis=0);den=np.sum(dr*dr,axis=1)
    f=np.clip(-np.sum(r[:-1]*dr,axis=1)/np.maximum(den,1e-300),0,1)
    d=np.linalg.norm(r[:-1]+f[:,None]*dr,axis=1);k=int(np.argmin(d))
    return {'minimum_mm':float(d[k]),'time_ps':float(t[k]+f[k]*(t[k+1]-t[k])),
            'meaning':'PIECEWISE_LINEAR_ESTIMATE_NOT_AN_ERROR_BOUND'}


def audit(path):
    result={'scope':'FIRST_RETURN_WINDOW_NUMERICAL_GEOMETRY',
            'archive_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'physical_annihilation_probability':None,'cases':[]}
    with tarfile.open(path) as tar:
        manifest=json.load(tar.extractfile('manifest.json'))
        result['source_commit']=manifest['observed_source_commit']
        for case in manifest['cases']:
            with np.load(io.BytesIO(tar.extractfile(case['name']+'.npz').read()),allow_pickle=False) as a:
                pair=[]
                for q in (0,1):
                    mask=a['q']==q
                    t=a['t'][mask]*1000;x=np.column_stack([a[k][mask] for k in ('x','y','z')])
                    if np.any(np.diff(t)<0):raise ValueError('time-reversed samples')
                    keep=np.r_[np.diff(t)>0,True];pair.append((t[keep],x[keep]))
                bound=speed_cone_window(*pair[0],*pair[1])
                radius=case['close_return']['match_radius_mm']
                lower=bound['minimum_separation_lower_bound_mm'];upper=bound['minimum_separation_upper_bound_mm']
                status=('OUTSIDE_MATCHING_RADIUS_OVER_WINDOW' if lower>radius else
                        'INSIDE_MATCHING_RADIUS_AT_SOME_TIME' if upper<radius else 'UNRESOLVED')
                result['cases'].append({'name':case['name'],'status':status,'match_radius_mm':radius,
                                        'speed_cone':bound,'chord_estimate':chord_minimum(pair)})
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=audit(a.archive);a.output.write_text(json.dumps(r,indent=2)+'\n')
    for c in r['cases']:
        print(c['name'],c['status'],'chord_min_um',c['chord_estimate']['minimum_mm']*1000,
              'lower_bound_um',c['speed_cone']['minimum_separation_lower_bound_mm']*1000)
