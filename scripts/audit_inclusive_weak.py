"""Streaming inclusive-handoff audit and compact return windows; ROOT stays local.

Vertex conservation/handoff tests are distinct from piecewise-linear geometry.
No fitted string energy, full-worldsheet Noether closure or capture is claimed.
"""
import argparse,hashlib,json,math,tarfile
from pathlib import Path
import numpy as np
C=.299792458


def parse_vertices(lines):
    dec=[];hand=[]
    for line in lines:
        if 'QP_WEAK_JSON ' not in line:continue
        d=json.loads(line.split('QP_WEAK_JSON ',1)[1])
        if d['kind']=='decay':
            if d.get('schema')!=2:raise ValueError('expected inclusive schema 2')
            total=np.array(d['daughter_mev'])+np.sum([x['p4_mev'] for x in d['products']],axis=0)
            residual=float(np.max(abs(np.array(d['parent_mev'])-total)))
            if residual>1e-6 or abs(d['proper_ps']-d['threshold_ps'])>2e-9:raise ValueError('weak conservation or clock failed')
            d['independent_four_residual_mev']=residual;dec.append(d)
        if d['kind']=='handoff':
            if not d['history_preserved']:raise ValueError('handoff fingerprint changed')
            hand.append(d)
    if len(dec)!=len(hand):raise ValueError('missing/extra handoff')
    return dec,hand


def local_minima(pair):
    lo=max(t[0] for t,x in pair);hi=min(t[-1] for t,x in pair)
    grid=np.unique(np.concatenate([t for t,x in pair]));grid=grid[(grid>=lo)&(grid<=hi)]
    if len(grid)<3:raise ValueError('insufficient common return window')
    x=[np.column_stack([np.interp(grid,t,v[:,j]) for j in range(3)]) for t,v in pair];r=x[0]-x[1];dr=np.diff(r,axis=0);dt=np.diff(grid);den=(dr*dr).sum(axis=1)
    f=-(r[:-1]*dr).sum(axis=1)/np.maximum(den,1e-300);rows=[]
    for k in np.flatnonzero((f>0)&(f<1)):
        rows.append((float(grid[k]+f[k]*dt[k]),float(np.linalg.norm(r[k]+f[k]*dr[k]))))
    left=(r[1:-1]*dr[:-1]).sum(axis=1);right=(r[1:-1]*dr[1:]).sum(axis=1)
    norm2=(r*r).sum(axis=1)
    knots=(left<=0)&(right>=0)&((norm2[1:-1]<norm2[:-2])|(norm2[1:-1]<norm2[2:]))
    for k in np.flatnonzero(knots)+1:rows.append((float(grid[k]),float(np.linalg.norm(r[k]))))
    rows.sort();return [{'time_cm_ps':t,'separation_mm':d} for t,d in rows]


def root_window(path,case,dec,maximum_samples=200000):
    import uproot
    b=case['pair_boost_z'];g=1/math.sqrt(1-b*b)
    after=max(g*(d['time_ps']-b*d['position_mm'][2]/C) for d in dec)
    M=case['mass_gev'];p=M*case['beta_each_cm']/math.sqrt(1-case['beta_each_cm']**2)
    period=2*p*6.582119569e-13/(case['lambda_ev']*1e-9)**2
    lo=after+1e-7;hi=after+3*period;pieces=[];n=0;last={0:0.,1:0.}
    fields=['q','t','x','y','z','px','py','pz','mass','pdg','proper','charge']
    with uproot.open(path) as root:
        for a in root['steps'].iterate(fields,step_size='8 MB',library='np'):
            for q in (0,1):
                m=a['q']==q
                if m.any():last[q]=max(last[q],float(a['t'][m].max())*1000)
            tc=g*(a['t']*1000-b*a['z']/C);keep=(abs(a['pdg'])==9900017)&(tc>=lo-period*.05)&(tc<=hi)
            n+=int(keep.sum())
            if n>maximum_samples:raise ValueError('return-window sample budget exceeded; retain ROOT for dedicated audit')
            if keep.any():pieces.append({k:a[k][keep] for k in fields})
    if not pieces:raise ValueError('no neutral samples in requested window')
    a={k:np.concatenate([p[k] for p in pieces]) for k in fields};pair=[]
    for q in (0,1):
        m=(a['q']==q);t=g*(a['t'][m]*1000-b*a['z'][m]/C);x=np.column_stack([a[k][m] for k in ('x','y','z')]);x[:,2]=g*(x[:,2]-b*C*a['t'][m]*1000)
        if len(t)<3 or np.any(np.diff(t)<0):raise ValueError('insufficient or time-reversed neutral history')
        dedup=np.r_[np.diff(t)>0,True];pair.append((t[dedup],x[dedup]))
    minima=[v for v in local_minima(pair) if v['time_cm_ps']>lo]
    if not minima:raise ValueError('no bracketed post-neutral return; window is censored')
    return a,{'status':'RECORDED_FIRST_RETURN_ESTIMATE','first_return':minima[0],'returns_in_saved_window':minima,'last_endpoint_lab_times_ps':last,'sample_count':n,'meaning':'PIECEWISE_LINEAR_RECORDED_GEOMETRY_NOT_CAPTURE_OR_TRANSPORT_ERROR_BOUND'}


def audit_and_pack(out):
    state=json.loads((out/'manifest.json').read_text());rows=[];packed=[]
    for c in state['cases']:
        name=c['name'];r={'case':name,'channel':c['channel'],'field_tesla':c['field_tesla'],'mode':c['mode'],'level':c['level']}
        try:
            if c.get('execution',{}).get('returncode')!=0 or not (out/(name+'.root')).is_file():raise ValueError('transport incomplete')
            if hashlib.sha256((out/(name+'.mac')).read_bytes()).hexdigest()!=c['macro_sha256']:raise ValueError('macro hash mismatch')
            def checked_lines(log):
                for line in log:
                    if 'QuirkMissingPartner' in line or 'WeakOrphan' in line:raise ValueError('orphan endpoint was reported')
                    yield line
            with (out/(name+'.log')).open(errors='replace') as log:dec,hand=parse_vertices(checked_lines(log))
            if c['weak_enabled']:
                if len(dec)!=2 or {d['parent_pdg'] for d in dec}!={17,-17}:raise ValueError('expected two forced endpoint decays')
                if any(d['channel']!=c['channel'] for d in dec):raise ValueError('forced decay channel mismatch')
                for d in dec:
                    expected=c['plan'][0 if d['parent_pdg']>0 else 1]
                    if abs(d['threshold_ps']-expected['proper_ps'])>1e-12:raise ValueError('clock differs from input plan')
                    if abs(d['daughter_mass_mev']-c['neutral_mass_gev']*1000)>1e-6:raise ValueError('neutral mass differs from manifest')
                if hashlib.sha256((out/(name+'.weakplan')).read_bytes()).hexdigest()!=c['general_plan_sha256']:raise ValueError('weak plan hash mismatch')
                a,geometry=root_window(out/(name+'.root'),c,dec);dest=out/(name+'.npz');np.savez_compressed(dest,**a);packed.append(dest);r['geometry']=geometry
            elif dec or hand:raise ValueError('stable case unexpectedly decayed')
            r.update(status='PASS_VERTEX_AND_HANDOFF_CHECKS',decays=dec,handoffs=hand)
        except (ValueError,KeyError,OSError) as e:r.update(status='FAIL_OR_CENSORED',error=str(e))
        rows.append(r);print('Audited',name,r['status'],flush=True)
    refinements=[]
    for r in rows:
        if r['level']!=7 or 'geometry' not in r:continue
        other=next(x for x in rows if x['case']==r['case'].replace('level7','level8'))
        if 'geometry' not in other:continue
        a=r['geometry']['first_return'];b=other['geometry']['first_return'];refinements.append({'case':r['case'],'first_return_time_shift_ps':b['time_cm_ps']-a['time_cm_ps'],'first_miss_shift_nm':(b['separation_mm']-a['separation_mm'])*1e6})
    report={'scope':state['scope'],'source_commit':state.get('observed_source_commit'),'cases':rows,'refinements':refinements,'all_vertex_handoff_checks_pass':all(r['status']=='PASS_VERTEX_AND_HANDOFF_CHECKS' for r in rows),'physical_survival_probability':None,'full_string_Noether_validation':'NOT_ESTABLISHED_BY_THIS_AUDIT'}
    audit=out/'inclusive_audit.json';audit.write_text(json.dumps(report,indent=2)+'\n');target=out.with_name(out.name+'_compact.tar.gz')
    with tarfile.open(target,'w:gz') as tar:
        tar.add(audit,arcname=audit.name);tar.add(out/'manifest.json',arcname='manifest.json')
        for c in state['cases']:
            for suffix in ('.mac','.card','.weakplan'):
                p=out/(c['name']+suffix)
                if p.exists():tar.add(p,arcname=p.name)
        for p in packed:tar.add(p,arcname=p.name)
        for name in ('build.log','patch_check.log','patch_apply.log'):
            if (out/name).is_file():tar.add(out/name,arcname=name)
    if not report['all_vertex_handoff_checks_pass']:print('Some checks failed/censored: return the compact archive and retain ROOT/logs.',flush=True)
    return target

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output_dir',type=Path);a=p.parse_args();print(audit_and_pack(a.output_dir.resolve()))
