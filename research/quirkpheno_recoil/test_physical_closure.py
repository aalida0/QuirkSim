import math
import unittest
from scipy.integrate import quad
from quirksim.physical_closure import capture_annihilation_interval,centrifugal_action_floor


class PhysicalClosureTest(unittest.TestCase):
    def test_missing_triplet_sector_is_not_assigned_zero_annihilation(self):
        r=capture_annihilation_interval(.8,1.,1e-5,.25)
        self.assertAlmostEqual(r['lower'],.2);self.assertAlmostEqual(r['upper'],.8)
        self.assertEqual(capture_annihilation_interval(.8,1.,1e-5),{'lower':0.,'upper':.8,'spin_status':'UNSPECIFIED'})

    def test_no_capture_and_zero_age(self):
        self.assertEqual(capture_annihilation_interval(0.,1.,1e-5)['upper'],0.)
        r=capture_annihilation_interval(1.,0.,1e-5,1.)
        self.assertEqual(r['lower'],0.);self.assertEqual(r['upper'],0.)
        self.assertEqual(capture_annihilation_interval(1.,0.,1e-5)['upper'],0.)
        with self.assertRaises(ValueError):capture_annihilation_interval(1.,1.,-1.)

    def test_action_bound_against_independent_quadrature(self):
        m=225.;e=40.;ell=100.;k=.2;H=2*m+e;sigma=(774e-9)**2
        actual=quad(lambda r:math.sqrt(m*m+(ell/r)**2-((H+ k/r-sigma*r)/2)**2),1/m,ell/H)[0]
        bound=centrifugal_action_floor(m,e,ell,k)
        self.assertGreater(actual,bound['action_lower_bound'])
        self.assertLess(bound['log10_WKB_exponential_ceiling'],0.)
        self.assertEqual(centrifugal_action_floor(m,e,0.)['status'],'NOT_IN_LARGE_L_DOMAIN')
