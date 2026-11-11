"""Compare endpoint fits with classical Noether targets from weak vertices.

Zero B gives a full pair-plus-string target. Uniform Bz gives only axial J.
Residuals are missing field contributions OR numerical error, not measured string J.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.audit_quirkpheno_return_window import chord_minimum

HBARC=1.973269804e-13
C=.299792458
QB=.000299792458  # GeV/mm per tesla, for charge +e


def boost_z(v,b):
    g=1/np.sqrt(1-b*b);v=np.array(v,float);out=v.copy()
    out[0]=g*(v[0]-b*v[3]);out[3]=g*(v[3]-b*v[0]);return out


def rest_angular(tensor,momentum):
    beta=momentum[1:]/momentum[0];b2=beta@beta
    transform=np.eye(4)
    if b2:
        g=1/np.sqrt(1-b2);transform[0,0]=g
        transform[0,1:]=transform[1:,0]=-g*beta
        transform[1:,1:]+=(g-1)*np.outer(beta,beta)/b2
    rest=transform@tensor@transform.T
    return np.array([rest[2,3],rest[3,1],rest[1,2]])/HBARC


def analyze(archive):
    rows=[]
    with tarfile.open(archive) as tar:
        manifest=json.load(tar.extractfile('manifest.json'))
        for case in manifest['cases']:
            if case['window_fit_status']!='BOTH_NEUTRAL_AT_FIRST_RETURN':continue
            b=case['pair_boost_z'];pair=[];mom=[]
            with np.load(io.BytesIO(tar.extractfile(case['name']+'.npz').read()),allow_pickle=False) as a:
                for q in (0,1):
                    mask=a['q']==q
                    xyz=np.array([a[k][mask] for k in ('x','y','z')]).T
                    tt=a['t'][mask]*1000
                    events=np.array([boost_z([C*t,*x],b) for t,x in zip(tt,xyz)])
                    pp=np.array([a[k][mask]/1000 for k in ('px','py','pz')]).T
                    ee=np.sqrt((a['mass'][mask]/1000)**2+(pp*pp).sum(axis=1))
                    vectors=np.array([boost_z([e,*p],b) for e,p in zip(ee,pp)])
                    times=events[:,0]/C
                    if np.any(np.diff(times)<0):raise ValueError('nonmonotone CM history')
                    keep=np.r_[np.diff(times)>0,True]
                    pair.append((times[keep],events[keep,1:]));mom.append(vectors[keep])
            minimum=chord_minimum(pair);time=minimum['time_ps']
            x=np.array([[np.interp(time,v[0],v[1][:,j]) for j in range(3)] for v in pair])
            p=np.array([[np.interp(time,pair[i][0],mom[i][:,j]) for j in range(4)] for i in (0,1)])
            endpoint=np.cross(x,p[:,1:]).sum(axis=0)/HBARC
            relative=np.cross(x[0]-x[1],(p[0,1:]-p[1,1:])/2)/HBARC
            extrinsic=np.cross((x[0]+x[1])/2,p[:,1:].sum(axis=0))/HBARC
            assert np.linalg.norm(endpoint-relative-extrinsic)<.1
            tensor=np.zeros((4,4));pions=np.zeros(4);torque_z=0.;decays=0
            for line in tar.extractfile(case['name']+'.log').read().decode().splitlines():
                if 'QP_WEAK_JSON ' not in line:continue
                d=json.loads(line.split('QP_WEAK_JSON ',1)[1])
                if d['kind']!='decay':continue
                event=boost_z([C*d['time_ps'],*d['position_mm']],b)
                if event[0]/C>time:continue
                pion=boost_z(np.array(d['pion_mev'])/1000,b)
                tensor-=np.outer(event,pion)-np.outer(pion,event);pions+=pion
                charge=1 if d['parent_pdg']>0 else -1
                torque_z-=QB*case['field_tesla']*charge*(event[1]**2+event[2]**2)/2/HBARC
                decays+=1
            if decays!=2:raise ValueError('requires both prior weak vertices')
            pion_target=np.array([tensor[2,3],tensor[3,1],tensor[1,2]])/HBARC
            axial=torque_z+pion_target[2]
            row={'case':case['name'],'point':case['point'],'mode':case['mode'],'seed':case['seed'],'level':case['level'],
                'relative_L_vector_hbar':relative.tolist(),'endpoint_J_vector_hbar':endpoint.tolist(),
                'endpoint_extrinsic_vector_hbar':extrinsic.tolist(),
                'axial_total_target_hbar':axial,'axial_required_field_or_numerical_residual_hbar':axial-endpoint[2],
                'axial_external_impulse_hbar':torque_z}
            if case['field_tesla']==0:
                initial=np.array([2*manifest['mass_gev']/np.sqrt(1-case['beta_each_cm']**2),0,0,0])
                surviving=initial-pions
                total_rest=rest_angular(tensor,surviving)
                row.update(total_J_target_vector_hbar=pion_target.tolist(),
                    total_rest_J_target_hbar=total_rest.tolist(),
                    total_rest_J_norm_hbar=float(np.linalg.norm(total_rest)),
                    required_field_or_numerical_residual_hbar=(pion_target-endpoint).tolist(),
                    relative_residual_norm=float(np.linalg.norm(pion_target-endpoint)/max(1.,np.linalg.norm(pion_target))))
            rows.append(row)
    return {'scope':'CLASSICAL_CONSERVATION_TARGET_NOT_DIRECT_STRING_QUADRATURE',
        'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'cases':rows,
        'assumptions':['Fixed production CM, boost along uniform Bz, zero initial orbital angular momentum.',
            'Only classical Lorentz torque and logged weak emissions; no material or radiation.',
            'Zero-B target is full pair-plus-string angular momentum; 2T target is its axial component only.',
            'Residual may be string angular momentum or numerical nonconservation; fields were not saved.',
            'Spin and finite initialization seed corrections are not resolved at O(1) hbar.',
            'Orbital target is not a quantum partial-wave occupation.']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=analyze(a.archive);a.output.write_text(json.dumps(r,indent=2)+'\n');print('checked',len(r['cases']))
