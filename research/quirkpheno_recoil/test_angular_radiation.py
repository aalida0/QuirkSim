import unittest
from quirksim.angular_radiation import reachability


class AngularBudgetTest(unittest.TestCase):
    def test_triangle_limit_and_threshold(self):
        self.assertFalse(reachability(12.,10.,1.)['reachable_by_any_bounded_step_sequence'])
        self.assertTrue(reachability(11.,10.,1.)['reachable_by_any_bounded_step_sequence'])

    def test_exhaustive_unbiased_walk_obeys_hitting_bound(self):
        import itertools
        n=10;start=5;hits=0
        for steps in itertools.product((-1,1),repeat=n):
            position=start;hit=False
            for step in steps:
                position+=step;hit|=position<=1
            hits+=hit
        bound=reachability(start,n,1.)['unbiased_vector_walk_low_L_hitting_probability_upper_bound']
        self.assertLessEqual(hits/2**n,bound)
        self.assertLess(bound,1.)

    def test_no_energy_already_low_and_invalid(self):
        self.assertEqual(reachability(0.,0.,1.)['unbiased_vector_walk_low_L_hitting_probability_upper_bound'],1.)
        self.assertEqual(reachability(2.,0.,1.)['unbiased_vector_walk_low_L_hitting_probability_upper_bound'],0.)
        with self.assertRaises(ValueError):reachability(1.,1.,0.)
