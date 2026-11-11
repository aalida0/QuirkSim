"""Synthetic 0/1/2 decay and Python/C++ plan roundtrip tests, not transport evidence."""
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import numpy as np
import awkward as ak
import uproot
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT))
from af_quirkpheno_weak_pilot import prepare
from audit_quirkpheno_sampled_pion import audit
from quirksim.decays import sample_decay
class FixedN:
    def __init__(self,n):self.n=n;self.count=0
    def uniform(self,a,b):
        self.count+=1
        return self.n[2] if self.count==1 else math.atan2(self.n[1],self.n[0])%(2*math.pi)
with tempfile.TemporaryDirectory(prefix='qp-sampled-check-') as tmp:
    tmp=Path(tmp);out=tmp/'synthetic';state=prepare(out,sampled_pion_seeds=8)
    assert len(state['cases'])==20 and state['sampled_seeds']==list(range(1,9))
    exe=tmp/'parse'
    subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic',str(HERE/'native_plan_check.cc'),'-o',str(exe)],check=True)
    for case in state['cases']:
        case['execution']={'returncode':0,'timeout':False};case['root_exists']=True
        if case['weak_enabled']:
            actual=json.loads(subprocess.check_output([str(exe),case['plan_file']],text=True))
            assert actual==[[p['proper_ps'],*p['direction']] for p in case['plan']]
        columns={k:[] for k in ['q','t','proper','pdg','mass','charge','px','py','pz']};counts=[];records=[]
        for q in (0,1):
            mass=state['mass_gev'];p=np.array([98.2*(-1)**q,0,63.4]);E=math.hypot(mass,np.linalg.norm(p))
            planned=case['plan'][q] if case['weak_enabled'] else None
            td=planned['proper_ps']*E/mass if planned else math.inf
            decay=td<=131
            if decay:
                k=sample_decay('pion',p,mass,225.,FixedN(planned['direction']))
                d=k['daughter'];pi=k['pion']
                records += [{'kind':'decay','parent_track':q+1,'parent_pdg':17*(-1)**q,'daughter_pdg':9900017*(-1)**q,
                             'pion_pdg':211*(-1)**q,'time_ps':td,'proper_ps':planned['proper_ps'],'threshold_ps':planned['proper_ps'],
                             'parent_mev':[1000*E,*(1000*p)],'daughter_mev':(1000*d).tolist(),'pion_mev':(1000*pi).tolist(),
                             'parent_mass_mev':mass*1000,'daughter_mass_mev':225000.,'pion_mass_mev':139.57039},
                            {'kind':'handoff','parent_track':q+1,'daughter_pdg':9900017*(-1)**q,'history_preserved':True}]
            times=sorted(set([0,20,40,60,80,100,131]+([td] if decay else [])))
            counts.append(len(times))
            for t in times:
                after=decay and t>td
                momentum=d[1:] if after else p
                proper=planned['proper_ps']+(t-td)*225./d[0] if after else t*mass/E
                row={'q':q,'t':t/1000,'proper':proper/1000,'pdg':(9900017 if after else 17)*(-1)**q,
                     'mass':225000. if after else mass*1000,'charge':0 if after else (-1)**q,
                     'px':momentum[0]*1000,'py':momentum[1]*1000,'pz':momentum[2]*1000}
                for key,value in row.items():columns[key].append(value)
        with uproot.recreate(out/(case['name']+'.root')) as f:
            f['events']={'q_reason':ak.Array([[6,6]]),'q_te':ak.Array([[.131,.131]]),'q_steps':ak.Array([counts])}
            f['steps']={k:np.array(v) for k,v in columns.items()}
        (out/(case['name']+'.log')).write_text('\n'.join('QP_WEAK_JSON '+json.dumps(r) for r in records))
    state['observed_source_commit']=state['source_commit'];(out/'manifest.json').write_text(json.dumps(state))
    def pack(path):
        with tarfile.open(path,'w:gz') as archive:
            for p in out.iterdir():archive.add(p,arcname=out.name+'/'+p.name)
    path=tmp/'fixture.tar.gz';pack(path);r=audit(path)
    assert r['local_checks_pass'],r
    direct=audit(out)
    assert direct['local_checks_pass'] and direct['sampled_seed_decay_count_histogram']==r['sampled_seed_decay_count_histogram']
    assert all(r['sampled_seed_decay_count_histogram'][i]>0 for i in (0,1,2))
    for seed in range(1,9):
        a=next(c for c in state['cases'] if c['name']==f'sampled_seed{seed}_level8')
        b=next(c for c in state['cases'] if c['name']==f'antipodal_seed{seed}_level8')
        for p,q in zip(a['plan'],b['plan']):assert p['proper_ps']==q['proper_ps'] and np.allclose(p['direction'],-np.array(q['direction']))
    log=out/'sampled_seed1_level8.log';log.write_text(log.read_text().replace('"history_preserved": true','"history_preserved": false',1))
    pack(path);assert not audit(path)['local_checks_pass']
print('PASS sampled plans, exact C++ roundtrip, paired clocks/antipodes, synthetic 0/1/2 decays and corrupted handoff rejection')
