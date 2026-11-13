"""Streaming classical causal-contour conservation audit; no capture inference."""
import argparse,hashlib,json,tarfile
from collections import deque
from pathlib import Path
import numpy as np
C=.299792458
# Numerical budgets in native MeV/mm units, not O(hbar) quantum certification.
BUDGET={'energy_mev':1e-3,'momentum_mev':1e-3,'orbital_mev_mm':1e-3,'contour_mm':1e-8}


def audit_stream(path,field,weak,duration,expected_vertices=None):
    maxima={k:0. for k in BUDGET};n=0;vertices=set();reasons=set();pairs=[];before=None;previous_step=-1
    last=None;first_reference=None
    with path.open() as stream:
      for line in stream:
        d=json.loads(line);n+=1
        if d['schema']!=1:raise ValueError('unsupported Noether schema')
        def array(k,size):
            a=np.asarray(d[k],float)
            if a.shape!=(size,) or not np.isfinite(a).all():raise ValueError('nonfinite or malformed '+k)
            return a
        p=array('mechanical_p_mev',4)+array('emitted_p_mev',4)-array('initial_p_mev',4)
        j=array('mechanical_j_mev_mm',16)+array('emitted_j_mev_mm',16)-array('initial_j_mev_mm',16)
        if not np.isfinite(d['external_torque_z_mev_mm']):raise ValueError('nonfinite torque')
        if field==0 and d['external_torque_z_mev_mm']!=0:raise ValueError('external torque in zero field')
        for key in ('mechanical_j_mev_mm','emitted_j_mev_mm','initial_j_mev_mm'):
            tensor=array(key,16).reshape(4,4)
            if abs(tensor+tensor.T).max()>1e-8:raise ValueError('angular tensor is not antisymmetric')
        j[6]-=d['external_torque_z_mev_mm']  # flattened M^12 = Jz
        if first_reference is None:first_reference=(array('initial_p_mev',4),array('initial_j_mev_mm',16))
        elif not (np.array_equal(array('initial_p_mev',4),first_reference[0]) and np.array_equal(array('initial_j_mev_mm',16),first_reference[1])):raise ValueError('production reference changed')
        if d['step']<previous_step:raise ValueError('step ordering reversed')
        previous_step=d['step'];vertices.add(d['vertices']);reasons.add(d['reason'])
        maxima['energy_mev']=max(maxima['energy_mev'],abs(float(p[0])))
        if field==0:maxima['momentum_mev']=max(maxima['momentum_mev'],float(abs(p[1:]).max()))
        # All six orbital Lorentz generators at B=0; only Jz at uniform Bz.
        maxima['orbital_mev_mm']=max(maxima['orbital_mev_mm'],float(abs(j).max()) if field==0 else abs(float(j[6])))
        maxima['contour_mm']=max(maxima['contour_mm'],float(abs(array('contour_closure_mm',4)).max()))
        corrected=np.r_[p,j]
        if d['reason']=='before_vertex':before=(d['vertices'],corrected)
        if d['reason']=='after_vertex':
            if before is None or d['vertices']!=before[0]+1:raise ValueError('unpaired vertex snapshot')
            pairs.append(float(abs(corrected-before[1]).max()));before=None
        last=d
    expected=(2 if weak else 0) if expected_vertices is None else expected_vertices
    if expected not in (0,1,2):raise ValueError("invalid expected vertex count")
    if n<3 or last is None or last['reason']!='end_event':raise ValueError('missing/short/endless diagnostic stream')
    if last['vertices']!=expected or vertices!=set(range(expected+1)) or len(pairs)!=expected:raise ValueError('incomplete weak-vertex coverage')
    if before is not None:raise ValueError('incomplete vertex pair')
    times=np.asarray(last['endpoint_ct_mm'],float)/C
    if times.shape!=(2,) or not np.isfinite(times).all() or min(times)<duration:raise ValueError('premature endpoint termination')
    failed=[k for k,v in maxima.items() if v>BUDGET[k]]
    if pairs and max(pairs)>1e-6:failed.append('vertex_charge_jump')
    return dict(status='FAIL' if failed else 'PASS_SAMPLED_CLASSICAL_CHARGES',failed=failed,snapshots=n,max_residuals=maxima,vertex_pairs=len(pairs),max_vertex_jump_native=max(pairs,default=0.),end_lab_times_ps=times.tolist(),momentum_gate='ALL_COMPONENTS' if field==0 else 'NOT_CONSERVED_EXTERNAL_FIELD',angular_gate='ALL_SIX_LORENTZ_GENERATORS' if field==0 else 'AXIAL_JZ_WITH_EXTERNAL_TORQUE_ONLY')


def audit_and_pack(out):
    manifest=json.loads((out/'manifest.json').read_text());rows=[]
    for c in manifest['cases']:
        row={'case':c['name'],'field_tesla':c['field_tesla'],'level':c['level'],'channel':c['channel']}
        try:
            if c.get('execution',{}).get('returncode')!=0:raise ValueError('transport missing or failed')
            for suffix,key in (('.mac','macro_sha256'),('.weakplan','general_plan_sha256')):
                if key in c and hashlib.sha256((out/(c['name']+suffix)).read_bytes()).hexdigest()!=c[key]:raise ValueError('input hash mismatch')
            row.update(audit_stream(out/(c['name']+'.noether.jsonl'),c['field_tesla'],c['weak_enabled'],c['nominal_stop_ps']))
        except (ValueError,KeyError,OSError,TypeError) as e:row.update(status='UNKNOWN_INCOMPLETE',error=str(e))
        rows.append(row)
    pairs=[]
    for r in rows:
        if r['level']!=7 or 'max_residuals' not in r:continue
        fine=next(v for v in rows if v['case']==r['case'].replace('level7','level8'))
        if 'max_residuals' in fine:pairs.append({'coarse':r['case'],'fine':fine['case'],'coarse_max':r['max_residuals'],'fine_max':fine['max_residuals']})
    report=dict(scope='SAMPLED_CLASSICAL_CAUSAL_CONTOUR_CHARGES_NO_RADIATION_NO_CAPTURE_NO_QUANTUM_SPIN',budgets=BUDGET,cases=rows,refinements=pairs,all_pass=bool(rows) and all(r['status']=='PASS_SAMPLED_CLASSICAL_CHARGES' for r in rows),physical_survival_probability=None)
    p=out/'noether_audit.json';p.write_text(json.dumps(report,indent=2)+'\n')
    tails=[]
    if not report['all_pass']:
        for name in ['build.log','noether_patch_check.log']+[c['name']+'.log' for c in manifest['cases'] if c.get('execution',{}).get('returncode')!=0]:
            log=out/name
            if log.exists():
                with log.open(errors='replace') as f:tails.append(name+'\n'+''.join(deque(f,maxlen=60)))
        (out/'failure_tails.txt').write_text('\n'.join(tails))
    target=out.with_name(out.name+'_summary.tar.gz')
    with tarfile.open(target,'w:gz') as tar:
        for name in ('manifest.json','noether_audit.json','failure_tails.txt'):
            if (out/name).exists():tar.add(out/name,arcname=name)
        for c in manifest['cases']:
            for suffix in ('.mac','.card','.weakplan','.noether.jsonl'):
                p=out/(c['name']+suffix)
                if p.exists():tar.add(p,arcname=p.name)
    print('Conservation:',sum(r['status']=='PASS_SAMPLED_CLASSICAL_CHARGES' for r in rows),'/',len(rows),'Summary:',target,flush=True)
    return target


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output_dir',type=Path);a=p.parse_args();audit_and_pack(a.output_dir.resolve())
