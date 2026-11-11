"""Matched first-return geometry contrasts; no annihilation inference."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile
import numpy as np
from audit_quirkpheno_return_window import chord_minimum


def describe(values):
    v=np.asarray(values,float)
    return {'n':len(v),'median_nm':float(np.median(v)),
            'min_nm':float(v.min()),'max_nm':float(v.max())} if len(v) else None


def analyze(archive):
    rows=[]
    with tarfile.open(archive) as tar:
        manifest=json.load(tar.extractfile('manifest.json'))
        for c in manifest['cases']:
            pair=[];b=c['pair_boost_z'];g=1/np.sqrt(1-b*b)
            with np.load(io.BytesIO(tar.extractfile(c['name']+'.npz').read()),allow_pickle=False) as a:
                for q in (0,1):
                    mask=a['q']==q;t=a['t'][mask]*1000
                    x=np.column_stack([a[k][mask] for k in ('x','y','z')])
                    cm_t=g*(t-b*x[:,2]/.299792458)
                    x[:,2]=g*(x[:,2]-b*.299792458*t)
                    if np.any(np.diff(cm_t)<0):raise ValueError('nonmonotonic CM time')
                    keep=np.r_[np.diff(cm_t)>0,True];pair.append((cm_t[keep],x[keep]))
            fit=chord_minimum(pair);lo=max(p[0][0] for p in pair);hi=min(p[0][-1] for p in pair)
            if not lo<fit['time_ps']<hi:raise ValueError('edge-censored minimum')
            delta=abs(fit['minimum_mm']-c['return_minimum_cm_chord']['minimum_mm'])*1e6
            if delta>.1:raise ValueError('compact and full-search estimates differ by >0.1 nm')
            rows.append({k:c[k] for k in ('name','point','mode','seed','level','decays_by_first_return')})
            rows[-1].update(minimum_nm=fit['minimum_mm']*1e6,compact_vs_full_search_difference_nm=delta)
    fine={(r['point'],r['mode'],r['seed']):r for r in rows if r['level']==8}
    summaries=[]
    for point in manifest['points']:
        name=point['name'];stable=fine[name,'stable',0]['minimum_nm']
        neutral=[s for s in manifest['sampled_seeds'] if all(
            fine[name,mode,s]['decays_by_first_return']==2
            for mode in ('isotropic','collinear_forward','collinear_backward'))]
        contrasts=[]
        for mode in ('collinear_forward','collinear_backward'):
            iso=[fine[name,'isotropic',s]['minimum_nm'] for s in neutral]
            ctl=[fine[name,mode,s]['minimum_nm'] for s in neutral]
            contrasts.append({'control':mode,'matched_neutral_seeds':neutral,
                'isotropic':describe(iso),'collinear':describe(ctl),
                'paired_isotropic_minus_collinear':describe(np.array(iso)-ctl),
                'paired_collinear_minus_stable':describe(np.array(ctl)-stable)})
        categories=[]
        for mode in ('isotropic','collinear_forward','collinear_backward'):
            for n in (0,1,2):
                subset=[r for r in fine.values() if r['point']==name and r['mode']==mode and r['decays_by_first_return']==n]
                categories.append({'mode':mode,'decays':n,'seeds':[r['seed'] for r in subset],
                                   'minimum':describe([r['minimum_nm'] for r in subset])})
        summaries.append({'point':name,'stable_minimum_nm':stable,'neutral_contrasts':contrasts,'categories':categories})
    field=[]
    for mode in ('stable','isotropic','collinear_forward','collinear_backward'):
        seeds=[0] if mode=='stable' else manifest['sampled_seeds']
        for s in seeds:
            a=fine['benchmark',mode,s];b=fine['zero_field',mode,s]
            field.append({'mode':mode,'seed':s,'decays_2T':a['decays_by_first_return'],
                'decays_0T':b['decays_by_first_return'],
                'minimum_2T_nm':a['minimum_nm'],'minimum_0T_nm':b['minimum_nm'],
                'paired_2T_minus_0T_nm':a['minimum_nm']-b['minimum_nm']})
    shifts=[]
    for r in rows:
        if r['level']==7:
            shifts.append({'point':r['point'],'mode':r['mode'],'seed':r['seed'],
                'level8_minus_level7_nm':fine[r['point'],r['mode'],r['seed']]['minimum_nm']-r['minimum_nm']})
    return {'scope':'PAIRED_FIRST_RETURN_GEOMETRY_MECHANISM_COMPARISON',
        'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
        'points':summaries,'matched_field_contrasts':field,'refinement':shifts,'cases':rows,
        'limitations':['Chord minima are estimates, not certified trajectory bounds.',
            'Stable-to-collinear includes mass change, longitudinal recoil and charge history.',
            'No charge-only intervention exists in this archive.',
            'Collinear-to-isotropic changes conditional decay directions; distances are nonlinear.',
            '0/1 decay rows are geometry only; no mixed-charge capture model is applied.',
            'No inclusive annihilation probability or causal additive decomposition is inferred.']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.write_text(json.dumps(analyze(a.archive),indent=2)+'\n')
