import hashlib
import math
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from af_triplet_recoil_ensemble import prepare,first_minimum
from compare_repeated_encounters import paired_control_summary,fitted_angular_momentum


class EnsemblePreparationTest(unittest.TestCase):
    def test_only_symmetry_protected_control_projects_angular_floor(self):
        case={'field_tesla':0.,'direction_mode':'collinear_cm','theta_cm':math.pi/2}
        r=np.array([0.,0.,1e-12]);p=np.array([100.,0.,1e-10])
        ell,d=fitted_angular_momentum(case,r,p)
        self.assertEqual(ell,0.);self.assertGreater(d['raw_endpoint_L_over_hbar'],100.)
        ell,d=fitted_angular_momentum(dict(case,field_tesla=2.),r,p)
        self.assertGreater(ell,100.);self.assertFalse(d['symmetry_projection_applied'])
        with self.assertRaises(ValueError):fitted_angular_momentum(case,np.array([0.,1e-4,0.]),p)

    def test_paired_summary_keeps_only_matched_neutral_cases(self):
        fits=[{'case':name,'mode':mode,'point':'point','seed':seed,'level':8}
              for name,mode,seed in [('a','isotropic',1),('b','collinear_forward',1),('c','isotropic',2)]]
        runs=[]
        for name,prob in [('a',.2),('b',.1),('c',.8)]:
            runs.append({'case':name,'control':'fitted_recoil',
                'path':{'model':{'epsilon_per_full_period':.1,'angular_loss_per_effective_emission':0.},'horizon_ps':1000.},
                'closure_scenarios':[{'status':'DECLARED_MODEL_ONLY','singlet_fraction_at_capture':.25,
                    'excited_optical_depth_per_eligible_encounter':0.,'annihilation_probability':prob}]})
        summary=paired_control_summary(fits,runs)
        self.assertEqual(len(summary),1);self.assertEqual(summary[0]['matched_seeds'],[1])
        self.assertAlmostEqual(summary[0]['mean_isotropic_minus_control_probability'],.1)
        self.assertIsNone(summary[0]['paired_standard_error'])

    def test_grid_pairs_and_lorentz_kinematics(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=prepare(Path(tmp)/'ensemble',seeds=2)
            self.assertEqual(len(d['cases']),44)
            cases=d['cases']
            for c in cases:
                p=Path(c['macro']);self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),c['macro_sha256'])
                rows=[s.split()[1:4] for s in p.read_text().splitlines() if s.startswith('/qsim/gun/p')]
                momenta=np.array(rows,float);mass=d['mass_gev'];energy=np.sqrt(mass*mass+np.sum(momenta*momenta,axis=1))
                inv=energy.sum()**2-np.sum(momenta.sum(axis=0)**2)
                self.assertAlmostEqual(inv/(4*mass*mass/(1-c['beta_each_cm']**2)),1.,places=12)
                if c['weak_enabled']:
                    same=[x for x in cases if x['point']==c['point'] and x['seed']==c['seed'] and x['level']==c['level'] and x['weak_enabled']]
                    self.assertEqual(len(same),3)
                    self.assertTrue(all(x['plan']==c['plan'] for x in same))
                    self.assertTrue(all(hashlib.sha256(Path(x['plan_file']).read_bytes()).hexdigest()==x['plan_sha256'] for x in same))
            self.assertTrue(d['conditional_control_not_isotropic'])

    def test_first_return_search_and_edge_censor(self):
        t=np.linspace(160,280,1001);a={};a['q']=np.repeat([0,1],len(t));a['t']=np.tile(t/1000,2)
        a['x']=np.r_[.1*(t-223),-.1*(t-223)];a['y']=np.zeros(2*len(t));a['z']=np.tile(.25*.299792458*t,2)
        minimum,labs=first_minimum(a,{'predicted_first_return_ps':223.,'pair_boost_z':.25})
        self.assertAlmostEqual(minimum['minimum_mm'],0.,places=12)
        self.assertTrue(all(abs(x-223.)<1e-10 for x in labs))
        a['x']=np.r_[.1*(t-100),-.1*(t-100)]
        with self.assertRaises(ValueError):first_minimum(a,{'predicted_first_return_ps':223.,'pair_boost_z':.25})


if __name__=='__main__':unittest.main()
