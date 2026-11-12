"""Inventory late weak decays and export first post-neutral local-return windows.

Runs on AF ROOT outputs; only compact windows/logs are transferred. Censored
events stay in population accounting, never become zero-annihilation events.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import numpy as np


def postneutral_minima(a,case):
    b=case['pair_boost_z'];g=1/np.sqrt(1-b*b);pair=[]
    for q in (0,1):
        mask=(a['q']==q)&(abs(a['pdg'])==9900017)
        t=a['t'][mask]*1000;x=np.array([a[k][mask] for k in ('x','y','z')]).T
        if len(t)<2:return None,[]
        cm_t=g*(t-b*x[:,2]/.299792458);x[:,2]=g*(x[:,2]-b*.299792458*t)
        if np.any(np.diff(cm_t)<0):raise ValueError('time-reversed neutral history')
        keep=np.r_[np.diff(cm_t)>0,True];pair.append((cm_t[keep],x[keep]))
    lo=max(v[0][0] for v in pair);hi=min(v[0][-1] for v in pair)
    if hi<=lo:return None,[]
    t=np.unique(np.r_[pair[0][0],pair[1][0]]);t=t[(t>=lo)&(t<=hi)]
    x=[np.array([np.interp(t,v[0],v[1][:,j]) for j in range(3)]).T for v in pair]
    return [float(lo),float(hi)],grid_minima(t,x[0]-x[1])


def grid_minima(t,r):
    dr=np.diff(r,axis=0);den=(dr*dr).sum(axis=1)
    f=-(r[:-1]*dr).sum(axis=1)/np.maximum(den,1e-300);minima=[]
    for i in np.where((f>0)&(f<1)&(den>0))[0]:
        minima.append({'time_cm_ps':float(t[i]+f[i]*(t[i+1]-t[i])),
                       'separation_mm':float(np.linalg.norm(r[i]+f[i]*dr[i]))})
    for i in range(1,len(t)-1):
        if r[i]@dr[i-1]<=0 and r[i]@dr[i]>=0 and (np.linalg.norm(dr[i-1])>0 or np.linalg.norm(dr[i])>0):
            minima.append({'time_cm_ps':float(t[i]),'separation_mm':float(np.linalg.norm(r[i]))})
    minima.sort(key=lambda v:v['time_cm_ps'])
    return minima


def bound_index(a,t,right=False):
    lo=0;hi=len(a)
    while lo<hi:
        mid=(lo+hi)//2
        if a[mid,0]<t or right and a[mid,0]==t:lo=mid+1
        else:hi=mid
    return lo


def mapped_minima(pair,block_rows=4096):
    if any(len(v)<2 for v in pair):return None,[]
    start=max(v[0,0] for v in pair);end=min(v[-1,0] for v in pair)
    if end<=start:return None,[]
    lo=start;previous=None;minima=[]
    while lo<end:
        left=[max(0,bound_index(v,lo)-1) for v in pair]
        hi=min(end,*(v[min(len(v)-1,i+block_rows),0] for v,i in zip(pair,left)))
        if hi<=lo:raise ValueError('non-advancing neutral scan')
        chunks=[v[i:min(len(v),bound_index(v,hi,right=True)+1)] for v,i in zip(pair,left)]
        t=np.unique(np.r_[lo,hi,chunks[0][:,0],chunks[1][:,0]]);t=t[(t>=lo)&(t<=hi)]
        xyz=[np.column_stack([np.interp(t,v[:,0],v[:,j]) for j in (1,2,3)]) for v in chunks]
        r=xyz[0]-xyz[1];dr=np.diff(r,axis=0)
        if previous is not None and r[0]@previous<=0 and r[0]@dr[0]>=0 and (np.linalg.norm(previous)>0 or np.linalg.norm(dr[0])>0):
            minima.append({'time_cm_ps':float(lo),'separation_mm':float(np.linalg.norm(r[0]))})
        minima.extend(grid_minima(t,r));previous=dr[-1];lo=hi
    return [float(start),float(end)],minima


def streamed_window(tree,case,tempdir,step_size='8 MB',block_rows=4096):
    b=case['pair_boost_z'];g=1/np.sqrt(1-b*b)
    paths=[Path(tempdir)/f'neutral{q}.bin' for q in (0,1)];pending=[None,None]
    with paths[0].open('wb') as f0,paths[1].open('wb') as f1:
        for a in tree.iterate(['q','pdg','t','x','y','z'],step_size=step_size,library='np'):
            for q,f in enumerate((f0,f1)):
                ix=(a['q']==q)&(abs(a['pdg'])==9900017)
                if not np.any(ix):continue
                t=a['t'][ix]*1000;z=a['z'][ix]
                rows=np.column_stack([g*(t-b*z/.299792458),a['x'][ix],a['y'][ix],g*(z-b*.299792458*t)])
                if pending[q] is not None:rows=np.vstack([pending[q],rows])
                if not np.isfinite(rows).all() or np.any(np.diff(rows[:,0])<0):raise ValueError('invalid neutral coordinates')
                rows=rows[np.r_[np.diff(rows[:,0])>0,True]]
                rows[:-1].astype('<f8').tofile(f);pending[q]=rows[-1:]
        for q,f in enumerate((f0,f1)):
            if pending[q] is not None:pending[q].astype('<f8').tofile(f)
    if any(p.stat().st_size<64 for p in paths):return None,[],None
    pair=[np.memmap(p,dtype='<f8',mode='r',shape=(p.stat().st_size//32,4)) for p in paths]
    span,minima=mapped_minima(pair,block_rows);del pair
    if not minima:return span,minima,None
    target=minima[0]['time_cm_ps'];selected={}
    for a in tree.iterate(step_size=step_size,library='np'):
        cm=g*(a['t']*1000-b*a['z']/.299792458)
        keep=(abs(a['pdg'])==9900017)&(abs(cm-target)<=.5)
        if not np.any(keep):continue
        names=a.dtype.names if isinstance(a,np.ndarray) else a.keys()
        for k in names:selected.setdefault(k,[]).append(a[k][keep])
    return span,minima,{k:np.concatenate(v) for k,v in selected.items()}


def package(out):
    import uproot
    original=(out/'manifest.json').read_bytes();m=json.loads(original)
    if not (out/'local_audit.json').exists():
        from audit_quirkpheno_sampled_pion import audit as run_audit
        result=run_audit(out,progress=True)
        tmp=out/'local_audit.json.tmp';tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(out/'local_audit.json')
    audit=(out/'local_audit.json').read_bytes()
    if not json.loads(audit)['local_checks_pass']:raise ValueError('full ROOT/handoff audit must pass first')
    population={};inventory=[];cases=[];windows={}
    for case in m['cases']:
        if not case['weak_enabled']:continue
        print('Packing '+case['name'],flush=True)
        root=out/(case['name']+'.root')
        with tempfile.TemporaryDirectory(prefix='neutral_scan_',dir=out) as work:
            with uproot.open(root,object_cache=None,array_cache=None) as f:
                span,minima,a=streamed_window(f['steps'],case,work)
        status='QUALIFIED' if minima else 'RIGHT_CENSORED_NO_POSTNEUTRAL_MINIMUM'
        record={'case':case['name'],'point':case['point'],'mode':case['mode'],'seed':case['seed'],'level':case['level'],
            'status':status,'common_neutral_span_cm_ps':span,'local_minima':minima,
            'definition':'all interior local minima of equal-CM-time piecewise-linear separation; no trajectory-error certification'}
        inventory.append(record)
        if case['mode']=='isotropic' and case['level']==8:
            p=population.setdefault(case['point'],{'population_seeds':0,'postneutral_qualified':0,'right_censored':0})
            p['population_seeds']+=1;p['postneutral_qualified' if minima else 'right_censored']+=1
        if not minima:continue
        target=minima[0]['time_cm_ps'];b=case['pair_boost_z'];g=1/np.sqrt(1-b*b)
        cm=g*(a['t']*1000-b*a['z']/.299792458)
        keep=(abs(a['pdg'])==9900017)&(abs(cm-target)<=.5)
        spans=[]
        for q in (0,1):
            tt=cm[keep&(a['q']==q)];spans.append((tt.min(),tt.max()) if len(tt)>1 else (float('inf'),-float('inf')))
        if not max(s[0] for s in spans)<target<min(s[1] for s in spans):
            raise ValueError('neutral return window does not bracket the minimum: '+case['name'])
        names=a.dtype.names if isinstance(a,np.ndarray) else a.keys()
        buf=io.BytesIO();np.savez_compressed(buf,**{k:a[k][keep] for k in names})
        windows[case['name']+'.npz']=buf.getvalue();new=dict(case)
        new.update(window_fit_status='BOTH_NEUTRAL_AT_FIRST_POSTNEUTRAL_RETURN',selected_return=record['local_minima'][0])
        cases.append(new)
    m.update(suite='triplet_postneutral_returns',cases=cases,ensemble_summary=population,
        postneutral_packer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        original_manifest_sha256=hashlib.sha256(original).hexdigest(),postneutral_inventory=inventory,
        population_warning='qualified windows only in cases; every omitted history is explicitly censored in postneutral_inventory and population counts')
    target=out.parent/(out.name+'_postneutral.tar.gz')
    with tarfile.open(target,'w:gz') as tar:
        payload={'manifest.json':(json.dumps(m,indent=2)+'\n').encode(),'original_manifest.json':original,'local_audit.json':audit,**windows}
        for name,data in payload.items():
            info=tarfile.TarInfo(name);info.size=len(data);tar.addfile(info,io.BytesIO(data))
        for case in cases:
            for suffix in ('.log','.mac','.card','.plan'):
                f=out/(case['name']+suffix)
                if f.exists():tar.add(f,arcname=f.name)
    return target


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output_dir',type=Path);a=p.parse_args()
    print('Compact post-neutral archive:',package(a.output_dir.resolve()))
