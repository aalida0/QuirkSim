"""Local audit of the forced-pion lifetime pilot; not an inclusive physics result."""
import argparse
import hashlib
import io
import json
from pathlib import PurePosixPath, Path
from types import SimpleNamespace
import tarfile
import numpy as np
import uproot

class DirectoryArchive:
    """Read the same audit inputs directly on AF without packing ROOT files."""
    def __init__(self,path):self.path=Path(path).resolve();self.prefix=self.path.name+'/'
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def getnames(self):return [self.prefix+p.name for p in self.path.iterdir() if p.is_file()]
    def getmembers(self):return [SimpleNamespace(name=name) for name in self.getnames()]
    def extractfile(self,item):
        name=item if isinstance(item,str) else item.name
        if not name.startswith(self.prefix):raise ValueError('unexpected audit prefix')
        target=(self.path/name[len(self.prefix):]).resolve()
        if target.parent!=self.path:raise ValueError('audit path outside directory')
        return io.BytesIO(target.read_bytes())

def stream_endpoints(tree,step_size='8 MB',neutral_mass_gev=225.):
    states=[dict(rows=0,neutral=False,neutral_ok=True,finite=True,clock_order=True,time_order=True) for _ in (0,1)]
    branches=['q','t','proper','pdg','mass','charge','px','py','pz']
    for a in tree.iterate(branches,step_size=step_size,library='np'):
        for q,s in enumerate(states):
            mask=a['q']==q
            if not np.any(mask):continue
            times=a['t'][mask]*1000;proper=a['proper'][mask]*1000
            if s['rows']==0:s['initial_mass']=float(a['mass'][mask][0])
            else:
                s['clock_order'] &= bool(proper[0]>=s['last_proper']-2e-9)
                s['time_order'] &= bool(times[0]>=s['last_time'])
            s['clock_order'] &= bool(np.all(np.diff(proper)>=-2e-9))
            s['time_order'] &= bool(np.all(np.diff(times)>=0))
            s['finite'] &= all(bool(np.isfinite(a[k][mask]).all()) for k in ('t','proper','mass','px','py','pz'))
            s['finite'] &= bool(np.isfinite(times).all() and np.isfinite(proper).all())
            neutral=mask&(abs(a['pdg'])==9900017)
            s['neutral'] |= bool(np.any(neutral))
            s['neutral_ok'] &= bool(np.all(a['charge'][neutral]==0) and np.all(abs(a['mass'][neutral]-1000*neutral_mass_gev)<1e-6))
            s['rows']+=len(times);s['last_time']=float(times[-1]);s['last_proper']=float(proper[-1])
    return states


def audit(path,progress=False,step_size='8 MB'):
    cases={};histogram={0:0,1:0,2:0}
    is_directory=Path(path).is_dir()
    with (DirectoryArchive(path) if is_directory else tarfile.open(path)) as archive:
        member=next(m for m in archive.getmembers() if m.name.endswith('/manifest.json'))
        manifest=json.load(archive.extractfile(member));prefix=member.name.rsplit('/',1)[0]+'/'
        if manifest['suite'] not in ('sampled_pion','triplet_ensemble'):raise ValueError('This auditor requires sampled-pion or triplet-ensemble plans')
        for case in manifest['cases']:
            if progress:print('Auditing '+case['name'],flush=True)
            name=case['name'];checks={};out={'checks':checks,'decays':[],'endpoints':[]};cases[name]=out
            checks['execution']=case.get('execution',{}).get('returncode')==0 and not case.get('execution',{}).get('timeout',True)
            if prefix+name+'.root' not in archive.getnames():checks['root_present']=False;continue
            checks['root_present']=True
            checks['macro_hash']=hashlib.sha256(archive.extractfile(prefix+name+'.mac').read()).hexdigest()==case['macro_sha256']
            if case['weak_enabled']:
                plan=archive.extractfile(prefix+PurePosixPath(case['plan_file']).name).read()
                checks['plan_hash']=hashlib.sha256(plan).hexdigest()==case['plan_sha256']
            log=archive.extractfile(prefix+name+'.log').read().decode(errors='replace')
            checks['no_exception']='*** G4Exception' not in log
            records=[json.loads(l.split('QP_WEAK_JSON ',1)[1]) for l in log.splitlines() if 'QP_WEAK_JSON ' in l]
            decays=[r for r in records if r['kind']=='decay'];handoffs=[r for r in records if r['kind']=='handoff']
            out['decays']=decays;out['handoffs']=handoffs
            checks['handoff_count']=len(handoffs)==len(decays)
            checks['at_most_two_decays']=len(decays)<=2
            handle=Path(path)/(name+'.root') if is_directory else archive.extractfile(prefix+name+'.root')
            with uproot.open(handle,object_cache=None,array_cache=None) as root:
                e=root['events'].arrays(library='ak');summaries=stream_endpoints(root['steps'],step_size,manifest['neutral_mass_gev'])
            checks['one_event']=len(e)==1
            for q in (0,1):
                s=summaries[q];ds=[d for d in decays if (d['parent_pdg']>0)==(q==0)]
                hs=[h for h in handoffs if (h['daughter_pdg']>0)==(q==0)]
                did_decay=s['neutral'];count=int(did_decay)
                checks[f'q{q}_event_count']=len(ds)==count and len(hs)==count
                checks[f'q{q}_completion']=int(e['q_reason'][0][q])==6 and float(e['q_te'][0][q])*1000>=case.get('nominal_stop_ps',manifest['nominal_stop_ps'])
                checks[f'q{q}_rows']=s['rows']==int(e['q_steps'][0][q]) and s['rows']>0
                checks[f'q{q}_last_row_matches_event']=s['rows']>0 and abs(s['last_time']-float(e['q_te'][0][q])*1000)<1e-8
                checks[f'q{q}_finite']=s['finite'];checks[f'q{q}_clock_order']=s['clock_order'];checks[f'q{q}_time_order']=s['time_order']
                checks[f'q{q}_initial_mass']=abs(s.get('initial_mass',0)-manifest['mass_gev']*1000)<1e-6
                if did_decay:
                    checks[f'q{q}_neutral_mass_charge']=s['neutral_ok']
                if case['weak_enabled'] and not did_decay:
                    checks[f'q{q}_censoring']=s['rows']>0 and s['last_proper']<case['plan'][q]['proper_ps']+2e-9
                if not case['weak_enabled']:checks[f'q{q}_stable']=not did_decay
                out['endpoints'].append({'q':q,'decays':count,'last_time_ps':s.get('last_time'),'last_proper_ps':s.get('last_proper')})
                if len(ds)==1 and len(hs)==1:
                    d=ds[0];h=hs[0];planned=case['plan'][q]
                    checks[f'q{q}_threshold']=abs(d['threshold_ps']-planned['proper_ps'])<1e-10 and abs(d['proper_ps']-planned['proper_ps'])<=2e-9
                    checks[f'q{q}_identity']=d['parent_track']==h['parent_track'] and d['daughter_pdg']==h['daughter_pdg'] and d['daughter_pdg']==(9900017 if q==0 else -9900017) and d['pion_pdg']==(211 if q==0 else -211)
                    checks[f'q{q}_history']=h['history_preserved'] is True
                    p=np.array(d['parent_mev']);daughter=np.array(d['daughter_mev']);pion=np.array(d['pion_mev'])
                    checks[f'q{q}_conservation']=bool(np.max(abs(p-daughter-pion))<1e-6)
                    for label,v in [('parent',p),('daughter',daughter),('pion',pion)]:
                        checks[f'q{q}_{label}_shell']=bool(abs(v[0]*v[0]-v[1:]@v[1:]-d[label+'_mass_mev']**2)<.01)
                    beta=p[1:]/p[0];b2=float(beta@beta);gamma=p[0]/d['parent_mass_mev'];rest=pion[1:].copy()
                    if b2:rest+=((gamma-1)*float(beta@pion[1:])/b2-gamma*pion[0])*beta
                    expected=np.array(planned['direction'])
                    if case.get('direction_mode')=='collinear':
                        expected=(p[1:]/np.linalg.norm(p[1:]) if np.linalg.norm(p[1:]) else expected)*case['collinear_sign']
                    if case.get('direction_mode')=='collinear_cm':
                        bz=case['pair_boost_z'];gg=1/np.sqrt(1-bz*bz)
                        pcm=p.copy();picm=pion.copy()
                        pcm[0]=gg*(p[0]-bz*p[3]);pcm[3]=gg*(p[3]-bz*p[0])
                        picm[0]=gg*(pion[0]-bz*pion[3]);picm[3]=gg*(pion[3]-bz*pion[0])
                        bb=pcm[1:]/pcm[0];bb2=float(bb@bb);rest=picm[1:].copy()
                        if bb2:rest+=((pcm[0]/d['parent_mass_mev']-1)*float(bb@picm[1:])/bb2-pcm[0]/d['parent_mass_mev']*picm[0])*bb
                        expected=(pcm[1:]/np.linalg.norm(pcm[1:]) if np.linalg.norm(pcm[1:])>1e-14*d['parent_mass_mev'] else expected)*case['collinear_sign']
                    checks[f'q{q}_direction']=bool(np.linalg.norm(rest/np.linalg.norm(rest)-expected)<1e-8)
            out['local_checks_pass']=all(checks.values())
            if case.get('seed') and case['mode'].startswith('sampled_seed') and case['level']==8:
                count=sum(d['time_ps']<=case.get('nominal_stop_ps',manifest['nominal_stop_ps']) for d in decays)
                histogram[count]=histogram.get(count,0)+1
    return {'archive_sha256':None if is_directory else hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            'scope':manifest['scope'],'source_commit':manifest.get('observed_source_commit'),
            'adapter_sha256':manifest['adapter_sha256'],'cases':cases,
            'local_checks_pass':all(c.get('local_checks_pass',False) for c in cases.values()),
            'sampled_seed_decay_count_histogram':histogram,'pion_only_conditional':True,
            'trajectory_refinement':'PENDING_PAIRED_SEED1_TRACE_COMPARISON',
            'physical_survival_probability':None,
            'note':'Counts are diagnostics of this forced-pion fixed-kinematics finite-window pilot, not inclusive rates or survival.'}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive');p.add_argument('--output',required=True);a=p.parse_args()
    r=audit(a.archive,progress=Path(a.archive).is_dir())
    with open(a.output,'w') as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
    print('SAMPLED_PION_LOCAL_CHECKS_PASS' if r['local_checks_pass'] else 'SAMPLED_PION_LOCAL_CHECKS_FAIL')
