import unittest
import numpy as np
from quirksim.close_return import monitor, _bound


class CloseReturnTest(unittest.TestCase):
    def test_randomized_bound_on_asynchronous_linear_paths(self):
        rng=np.random.default_rng(14); n=10000
        a=rng.uniform(-1,1,(n,3)); b=rng.uniform(-1,1,(n,3))
        va=rng.uniform(-.05,.05,(n,3)); vb=rng.uniform(-.05,.05,(n,3))
        t0=rng.uniform(0,1,n); t1=t0+rng.uniform(0,1,n); tp=rng.uniform(-.5,1.5,n)
        ae=a+va*(t1-t0)[:,None]; bp=b+vb*tp[:,None]
        lower=_bound(a,t0,ae,t1,bp,tp)
        for f in np.linspace(0,1,51):
            t=t0+f*(t1-t0)
            distance=np.linalg.norm(a+va*(t-t0)[:,None]-b-vb*t[:,None],axis=1)
            self.assertTrue(np.all(lower<=distance+1e-12))

    def test_first_unproven_interval_contains_or_precedes_return(self):
        t=np.linspace(0.,50.,501)
        x0=np.column_stack((.1*t,np.zeros((len(t),2))))
        x1=np.column_stack((10.-.1*t,np.zeros((len(t),2))))
        r=monitor(t,x0,t,x1,start_ps=0,arm_radius_mm=1.,match_radius_mm=.1)
        self.assertEqual(r['status'],'CENSORED_UNRESOLVED_INTERVAL')
        self.assertLessEqual(r['censor_interval_ps'][0],49.5)
        self.assertGreater(r['armed_at_ps'], -1)

    def test_stable_separation_and_asynchronous_steps(self):
        t0=np.arange(0,10,.02);t1=np.arange(.005,10,.03)
        x0=np.column_stack((.1*t0,np.zeros((len(t0),2))))
        x1=np.column_stack((1+.1*t1,np.zeros((len(t1),2))))
        r=monitor(t0,x0,t1,x1,start_ps=.1,arm_radius_mm=.2,match_radius_mm=.01)
        self.assertEqual(r['status'],'OUTSIDE_RADIUS_UNTIL_COMMON_END')
        self.assertGreater(r['lower_bound_mm'],.9)

    def test_no_arming_and_discontinuous_handoff(self):
        t=[0,1,2]; x=np.zeros((3,3))
        self.assertEqual(monitor(t,x,t,x,start_ps=0,arm_radius_mm=1,match_radius_mm=.1)['status'],'NEVER_ARMED')
        with self.assertRaises(ValueError):
            monitor([0,1,1,2],[[0,0,0],[0,0,0],[1,0,0],[1,0,0]],t,x,
                    start_ps=0,arm_radius_mm=1,match_radius_mm=.1)

    def test_missing_coverage_and_bad_config(self):
        x=np.zeros((2,3))
        r=monitor([0,1],x,[2,3],x,start_ps=0,arm_radius_mm=1,match_radius_mm=.1)
        self.assertEqual(r['status'],'INSUFFICIENT_OVERLAPPING_HISTORY')
        with self.assertRaises(ValueError):
            monitor([0,1],x,[0,1],x,start_ps=0,arm_radius_mm=.1,match_radius_mm=1)


if __name__=='__main__': unittest.main()
