import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from package_triplet_postneutral import postneutral_minima,package
from complete_triplet_closure import cp_interval


class PostneutralTest(unittest.TestCase):
    def test_ROOT_compact_roundtrip_keeps_censored_population(self):
        import json,tarfile,tempfile
        import uproot
        t=np.linspace(0,100,1001);a={'q':np.repeat([0,1],len(t)),
            't':np.tile(t/1000,2),'pdg':np.tile(np.where(t>40,9900017,17),2),
            'x':np.r_[.01*(t-70),.01*(70-t)],'y':np.zeros(2*len(t)),'z':np.zeros(2*len(t))}
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'run';out.mkdir()
            cases=[{'name':name,'point':'benchmark','mode':'isotropic','seed':seed,'level':8,
                    'weak_enabled':True,'pair_boost_z':0.} for seed,name in [(1,'late'),(2,'censored')]]
            for c in cases:
                arr=dict(a)
                if c['name']=='censored':arr['pdg']=np.full(2*len(t),17)
                with uproot.recreate(out/(c['name']+'.root')) as root:root['steps']=arr
            (out/'manifest.json').write_text(json.dumps({'cases':cases}))
            (out/'local_audit.json').write_text(json.dumps({'local_checks_pass':True}))
            with tarfile.open(package(out)) as tar:
                manifest=json.load(tar.extractfile('manifest.json'))
                self.assertEqual(len(manifest['cases']),1)
                self.assertEqual(manifest['ensemble_summary']['benchmark'],
                    {'population_seeds':2,'postneutral_qualified':1,'right_censored':1})
                self.assertIn('late.npz',tar.getnames());self.assertNotIn('censored.npz',tar.getnames())

    def test_return_after_late_neutralization_and_censor(self):
        t=np.linspace(0,100,1001);pdg=np.where(t>40,9900017,17)
        a={'q':np.repeat([0,1],len(t)),'t':np.tile(t/1000,2),'pdg':np.tile(pdg,2),
            'x':np.r_[t-70,70-t],'y':np.zeros(2*len(t)),'z':np.zeros(2*len(t))}
        span,events=postneutral_minima(a,{'pair_boost_z':0.})
        self.assertEqual(len(events),1);self.assertAlmostEqual(events[0]['time_cm_ps'],70.)
        self.assertEqual(events[0]['separation_mm'],0.);self.assertGreater(span[0],40.)
        a['pdg']=np.full(2*len(t),17)
        self.assertEqual(postneutral_minima(a,{'pair_boost_z':0.}),(None,[]))

    def test_exact_binomial_edges_and_pilot(self):
        self.assertEqual(cp_interval(0,8)[0],0.);self.assertEqual(cp_interval(8,8)[1],1.)
        lo,hi=cp_interval(3,8);self.assertLess(lo,3/8);self.assertGreater(hi,3/8)
