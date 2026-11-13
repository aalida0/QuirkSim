"""Audit a censored inclusive LO weak ensemble; controls are paired, not events."""
import argparse,hashlib,json,math,tarfile
from collections import Counter,deque
from pathlib import Path
import numpy as np
from audit_inclusive_weak import parse_vertices
from audit_native_noether import audit_stream,BUDGET
C=.299792458
HC=1.973269804e-10  # MeV mm


def encounter_summary(path,both_neutral,after_cm):
    count=0;last_time=-math.inf;first=None;minimum=math.inf;end=None
    with path.open() as stream:
      for line in stream:
        d=json.loads(line)
        if end is not None:raise ValueError('encounter records after end marker')
        if d['kind']=='end_event':end=d;continue
        if d['kind']!='return':raise ValueError('unknown encounter record')
        if not both_neutral:raise ValueError('neutral return before two decays')
        t=float(d['time_cm_ps']);miss=float(d['separation_mm'])
        if not math.isfinite(t) or not math.isfinite(miss) or miss<0 or t<=last_time or t<after_cm-1e-9:raise ValueError('invalid return time or distance')
        last_time=t;points=d['endpoints']
        if len(points)!=2 or [v['pdg'] for v in points]!=[9900017,-9900017]:raise ValueError('return endpoints are not the neutral pair')
        x=np.array([v['x_mm'] for v in points],float);p=np.array([v['p4_mev'] for v in points],float)
        if x.shape!=(2,3) or p.shape!=(2,4) or not np.isfinite(x).all() or not np.isfinite(p).all():raise ValueError('malformed return phase space')
        if abs(np.linalg.norm(x[0]-x[1])-miss)>1e-9:raise ValueError('return distance inconsistent with positions')
        if first is None:
            first=dict(d,classical_relative_orbital_hbar=(np.cross(x[0]-x[1],(p[0,1:]-p[1,1:])/2)/HC).tolist(),classical_endpoint_J_hbar=(np.cross(x,p[:,1:]).sum(axis=0)/HC).tolist())
        minimum=min(minimum,miss);count+=1
    if end is None or end['returns']!=count or end['processed_intervals']<1:raise ValueError('incomplete encounter stream')
    return dict(status='RECORDED_POSTNEUTRAL_MINIMA' if count else ('CENSORED_NO_BRACKETED_RETURN' if both_neutral else 'CENSORED_NOT_BOTH_NEUTRAL'),count=count,first_return=first,minimum_saved_mm=None if not count else minimum,buffer_high_water=end['buffer_high_water'],meaning='COMMON_CM_PIECEWISE_LINEAR_CLASSICAL_GEOMETRY_NOT_CAPTURE_OR_A_QUANTUM_PARTIAL_WAVE')


def audit_case(out,c):
    name=c['name'];r={k:c[k] for k in ('name','seed','field_tesla','level','mode')}
    if c.get('execution',{}).get('returncode')!=0 or not c.get('root_exists'):raise ValueError('transport incomplete')
    if hashlib.sha256((out/(name+'.mac')).read_bytes()).hexdigest()!=c['macro_sha256']:raise ValueError('macro hash mismatch')
    def checked(log):
        for line in log:
            if 'QuirkMissingPartner' in line or 'WeakOrphan' in line:raise ValueError('orphan endpoint')
            yield line
    with (out/(name+'.log')).open(errors='replace') as log:dec,hand=parse_vertices(checked(log))
    if len(dec)>2 or len({v['parent_pdg'] for v in dec})!=len(dec):raise ValueError('duplicate/extra endpoint decay')
    with (out/(name+'.noether.jsonl')).open() as f:
        last=deque(f,maxlen=1)
    if not last:raise ValueError('empty Noether stream')
    end=json.loads(last[0]);proper=np.asarray(end['endpoint_proper_ps'],float)
    if end['reason']!='end_event' or proper.shape!=(2,) or not np.isfinite(proper).all():raise ValueError('missing terminal proper clocks')
    if c['weak_enabled']:
        if hashlib.sha256((out/(name+'.weakplan')).read_bytes()).hexdigest()!=c['general_plan_sha256']:raise ValueError('weak plan hash mismatch')
        for i,p in enumerate(c['plan']):
            selected=[v for v in dec if v['parent_pdg']==(17 if i==0 else -17)]
            expected_pdg=(9900017 if i==0 else -9900017) if selected else (17 if i==0 else -17)
            if end['endpoint_pdg'][i]!=expected_pdg:raise ValueError('terminal endpoint identity mismatch')
            if selected:
                v=selected[0]
                if v['channel']!=p['channel'] or abs(v['threshold_ps']-p['proper_ps'])>1e-10:raise ValueError('sampled input was not replayed')
                if proper[i]+2e-9<p['proper_ps']:raise ValueError('decay beyond terminal proper clock')
            elif p['proper_ps']<=proper[i]-2e-9:raise ValueError('decay threshold crossed but no decay was recorded')
    elif dec or hand:raise ValueError('stable control decayed')
    b=c['pair_boost_z'];g=1/math.sqrt(1-b*b)
    after=max((g*(v['time_ps']-b*v['position_mm'][2]/C) for v in dec),default=0.)
    geometry=encounter_summary(out/(name+'.encounters.jsonl'),len(dec)==2,after)
    conservation=audit_stream(out/(name+'.noether.jsonl'),c['field_tesla'],c['weak_enabled'],c['nominal_stop_ps'],expected_vertices=len(dec))
    r.update(status='PASS_CLASSICAL_PILOT' if conservation['status']=='PASS_SAMPLED_CLASSICAL_CHARGES' else 'CONSERVATION_BUDGET_FAILED',decay_count=len(dec),decays=dec,handoffs=hand,terminal_proper_ps=proper.tolist(),geometry=geometry,conservation=conservation)
    return r


def summarize(state,rows):
    cohorts=[];comparisons=[];refinements=[]
    for field in (0.,2.):
        selected=[r for r in rows if r['mode']=='recoil' and r['level']==8 and r['field_tesla']==field]
        counts=Counter(str(r.get('decay_count','UNKNOWN')) for r in selected)
        channels=Counter(v['channel'] for r in selected for v in r.get('decays',[]))
        cohorts.append(dict(field_tesla=field,independent_events=len(selected),qualified_events=sum(r['status']=='PASS_CLASSICAL_PILOT' for r in selected),decay_counts=dict(counts),observed_decay_channels=dict(channels),note='Finite-window, fixed-production pilot. Controls/refinements are not additional independent events.'))
        for seed in state['sampled_seeds']:
            paired=[next((r for r in rows if r['seed']==seed and r['field_tesla']==field and r['level']==8 and r['mode']==mode),None) for mode in ('recoil','collinear_forward','collinear_backward')]
            if any(r is None for r in paired):
                comparisons.append(dict(seed=seed,field_tesla=field,status='NOT_RUN_REQUIRED_CONTROLS'));continue
            r=dict(seed=seed,field_tesla=field,status='CENSORED_OR_UNQUALIFIED',case_statuses=[v['status'] for v in paired])
            if all(v['status']=='PASS_CLASSICAL_PILOT' and v['geometry']['first_return'] is not None for v in paired):
                v=[x['geometry']['first_return'] for x in paired]
                r.update(status='QUALIFIED_PAIRED_GEOMETRY',first_miss_nm={mode:x['separation_mm']*1e6 for mode,x in zip(('recoil','collinear_forward','collinear_backward'),v)},time_cm_ps={mode:x['time_cm_ps'] for mode,x in zip(('recoil','collinear_forward','collinear_backward'),v)},ordinal_matching='First minimum after both neutral; not necessarily the same oscillation index.')
            comparisons.append(r)
    for r in rows:
        if r['level']!=7:continue
        fine=next(v for v in rows if v['name']==r['name'].replace('level7','level8'))
        a=r.get('geometry',{}).get('first_return');b=fine.get('geometry',{}).get('first_return')
        if a and b:refinements.append(dict(coarse=r['name'],fine=fine['name'],first_miss_shift_nm=(b['separation_mm']-a['separation_mm'])*1e6,first_time_shift_ps=b['time_cm_ps']-a['time_cm_ps']))
    return dict(scope=state['scope'],cases=rows,cohorts=cohorts,paired_comparisons=comparisons,refinements=refinements,weak_reference=state['weak_reference'],conservation_budgets=BUDGET,all_transport_and_conservation_pass=bool(rows) and all(r['status']=='PASS_CLASSICAL_PILOT' for r in rows),physical_survival_probability=None,capture_model='NONE',annihilation_prediction='NOT_COMPUTED')


def audit_and_pack(out):
    state=json.loads((out/'manifest.json').read_text());rows=[]
    for c in state['cases']:
        try:r=audit_case(out,c)
        except (ValueError,KeyError,OSError,TypeError) as e:r={k:c[k] for k in ('name','seed','field_tesla','level','mode')};r.update(status='UNKNOWN_INCOMPLETE',error=str(e))
        rows.append(r);print('Audited',c['name'],r['status'],flush=True)
    report=summarize(state,rows);audit=out/'physical_weak_audit.json';audit.write_text(json.dumps(report,indent=2)+'\n')
    failed=[r['name']+'.log' for r in rows if r['status']!='PASS_CLASSICAL_PILOT'];tails=[]
    for name in ['build.log','noether_patch_check.log']+failed:
        p=out/name
        if p.exists():
            with p.open(errors='replace') as f:tails.append(name+'\n'+''.join(deque(f,maxlen=40)))
    tail=out/'failure_tails.txt';tail.write_text('\n'.join(tails))
    target=out.with_name(out.name+'_summary.tar.gz')
    with tarfile.open(target,'w:gz') as tar:
        for p in (out/'manifest.json',audit,tail):tar.add(p,arcname=p.name)
        for c in state['cases']:
            for suffix in ('.mac','.card','.weakplan','.encounters.jsonl'):
                p=out/(c['name']+suffix)
                if p.exists():tar.add(p,arcname=p.name)
    print('Compact ensemble:',target,'bytes:',target.stat().st_size,flush=True);return target


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output_dir',type=Path);a=p.parse_args();audit_and_pack(a.output_dir.resolve())
