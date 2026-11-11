import sys, pathlib,json,tarfile,tempfile,shutil
from check_native import check
import numpy as np,awkward as ak,uproot
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]/'scripts'))
from audit_quirkpheno_weak_pilot import audit
from af_quirkpheno_weak_pilot import prepare
temporary=tempfile.TemporaryDirectory(prefix='qp-audit-fixture-')
root=pathlib.Path(temporary.name)
out=root/'synthetic'
state=prepare(out);state['observed_source_commit']=state['source_commit']
fixtures=check()
for case in state['cases']:
 case['execution']={'returncode':0,'timeout':False};case['root_exists']=True
 records=[];columns={k:[] for k in ['q','t','x','y','z','px','py','pz','ekin','mass','charge','proper','pdg','track','parent','handoff']}
 counts=[];initials=[]
 for q in [0,1]:
  k=next(r for r in fixtures if r['sign']==case['pion_rest_y_sign'] and np.sign(r['parent'][1])==(-1)**q)
  p=np.array(k['parent']);d=np.array(k['daughter']);pi=np.array(k['pion']);initials.append(p);n=0
  decays=case['weak_enabled'] and q==case['endpoint'];td=60.;xd=.299792458*p[1:]/p[0]*td
  if decays:
   records=[{'kind':'decay','parent_track':q+1,'parent_pdg':17*(-1)**q,'daughter_pdg':9900017*(-1)**q,
     'pion_pdg':211*(-1)**q,'time_ps':td,'proper_ps':50.001,'threshold_ps':50.,'position_mm':xd.tolist(),
     'parent_mass_mev':225157.1702763212,'daughter_mass_mev':225000.,'pion_mass_mev':139.57039,
     'parent_mev':p.tolist(),'daughter_mev':d.tolist(),'pion_mev':pi.tolist()},
     {'kind':'handoff','parent_track':q+1,'daughter_track':3,'daughter_pdg':9900017*(-1)**q,'history_preserved':True}]
  for t in ([0,20,40,60,80,100,131] if not decays else [0,20,40,60,80,100,131]):
   after=decays and t>td;handoff=decays and t==td
   v=d if after else p;mass=225000. if after else 225157.1702763212
   x=(xd+.299792458*d[1:]/d[0]*(t-td)) if after else .299792458*p[1:]/p[0]*t
   row={'q':q,'t':t/1000,'x':x[0],'y':x[1],'z':x[2],
        'px':v[1],'py':v[2],'pz':v[3],'ekin':v[0]-mass,'mass':mass,
        'charge':0 if after else (-1)**q,'proper':t/1000,'pdg':(9900017 if after else 17)*(-1)**q,
        'track':3 if after else q+1,'parent':q+1 if after else 0,'handoff':int(handoff)}
   for k,v in row.items():columns[k].append(v)
   n+=1
  counts.append(n)
 events={'q_reason':ak.Array([[6,6]]),'q_te':ak.Array([[.131,.131]]),'q_steps':ak.Array([counts]),
         'q_e0':ak.Array([[x[0] for x in initials]])}
 for i,k in enumerate(['q_px0','q_py0','q_pz0'],1):events[k]=ak.Array([[x[i] for x in initials]])
 for k in ['q_x0','q_y0','q_z0']:events[k]=ak.Array([[0.,0.]])
 with uproot.recreate(out/(case['name']+'.root')) as f:f['events']=events;f['steps']={k:np.array(v) for k,v in columns.items()}
 (out/(case['name']+'.log')).write_text('\n'.join('QP_WEAK_JSON '+json.dumps(r) for r in records))
(out/'manifest.json').write_text(json.dumps(state))
def pack(path):
 with tarfile.open(path,'w:gz') as t:
  for p in out.iterdir():t.add(p,arcname=out.name+'/'+p.name)
pack(str(root/'good.tar.gz'));r=audit(str(root/'good.tar.gz'))
assert r['local_handoff_checks_pass'],r
assert len(r['resolution_comparisons'])==8 and len(r['paired_direction_comparisons'])==2
name='pion_yplus_end0_level7';p=out/(name+'.log');p.write_text(p.read_text().replace('"history_preserved": true','"history_preserved": false'))
pack(str(root/'bad.tar.gz'));bad=audit(str(root/'bad.tar.gz'));assert not bad['local_handoff_checks_pass']
print('PASS: synthetic eight-case ROOT/log audit; piecewise interpolation; corrupted handoff rejected (synthetic fixtures are not transport validation)')

temporary.cleanup()
