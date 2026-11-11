import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import awkward as ak
import uproot
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from af_triplet_recoil_ensemble import prepare
from audit_quirkpheno_sampled_pion import audit
from quirksim.decays import sample_decay,boost


class FixedN:
    def __init__(self,n):self.n=n;self.count=0
    def uniform(self,a,b):
        self.count+=1
        return self.n[2] if self.count==1 else math.atan2(self.n[1],self.n[0])%(2*math.pi)


class EnsembleAuditTest(unittest.TestCase):
    def test_exact_cm_collinear_records_and_corruption_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'data';state=prepare(out,1)
            c=next(c for c in state['cases'] if c['point']=='benchmark' and c['mode']=='collinear_backward' and c['level']==8)
            state['cases']=[c];state['observed_source_commit']=state['source_commit']
            c['execution']={'returncode':0,'timeout':False};M=state['mass_gev'];mass=225.;b=np.array([0,0,.25])
            columns={k:[] for k in ('q','t','proper','pdg','mass','charge','px','py','pz')};records=[];counts=[]
            for q in (0,1):
                pcm=np.array([(-1)**q*M*.4/math.sqrt(.84),0,0]);parent=boost(math.hypot(M,np.linalg.norm(pcm)),pcm,b)
                direction=-pcm/np.linalg.norm(pcm)
                k=sample_decay('pion',pcm,M,mass,FixedN(direction))
                daughter=boost(k['daughter'][0],k['daughter'][1:],b);pion=boost(k['pion'][0],k['pion'][1:],b)
                tau=c['plan'][q]['proper_ps'];td=tau*parent[0]/M
                records.extend([{'kind':'decay','parent_track':q+1,'parent_pdg':17*(-1)**q,'daughter_pdg':9900017*(-1)**q,'pion_pdg':211*(-1)**q,
                    'threshold_ps':tau,'proper_ps':tau,'time_ps':td,'parent_mev':(1000*parent).tolist(),'daughter_mev':(1000*daughter).tolist(),
                    'pion_mev':(1000*pion).tolist(),'parent_mass_mev':1000*M,'daughter_mass_mev':1000*mass,'pion_mass_mev':139.57039},
                    {'kind':'handoff','parent_track':q+1,'daughter_pdg':9900017*(-1)**q,'history_preserved':True}])
                times=sorted(set([0.,td,td+1.,c['nominal_stop_ps']+1.]));counts.append(len(times))
                for t in times:
                    after=t>td;p=daughter if after else parent
                    row={'q':q,'t':t/1000,'proper':(tau+(t-td)*mass/daughter[0] if after else t*M/parent[0])/1000,
                         'pdg':(9900017 if after else 17)*(-1)**q,'mass':1000*(mass if after else M),'charge':0 if after else (-1)**q,
                         'px':1000*p[1],'py':1000*p[2],'pz':1000*p[3]}
                    for key,value in row.items():columns[key].append(value)
            with uproot.recreate(out/(c['name']+'.root')) as root:
                root['events']={'q_steps':ak.Array([counts]),'q_reason':ak.Array([[6,6]]),'q_te':ak.Array([[(c['nominal_stop_ps']+1)/1000]*2])}
                root['steps']={k:np.array(v) for k,v in columns.items()}
            (out/(c['name']+'.log')).write_text('\n'.join('QP_WEAK_JSON '+json.dumps(r) for r in records))
            (out/'manifest.json').write_text(json.dumps(state))
            result=audit(out);self.assertTrue(result['local_checks_pass'],result)
            records[0]['pion_mev'][1]+=.01
            (out/(c['name']+'.log')).write_text('\n'.join('QP_WEAK_JSON '+json.dumps(r) for r in records))
            self.assertFalse(audit(out)['local_checks_pass'])


if __name__=='__main__':unittest.main()
