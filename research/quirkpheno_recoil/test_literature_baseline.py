import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from literature_deexcitation_baseline import eligibility
from quirksim.radiation_annihilation import OrbitAveragedModel,linear_cooling_time_ps


class LiteratureBaselineTest(unittest.TestCase):
    def test_three_valued_gate_and_invalid_bounds(self):
        self.assertEqual(eligibility(2.,3.,1.),[0.,0.])
        self.assertEqual(eligibility(.1,.2,1.),[1.,1.])
        self.assertEqual(eligibility(0.,2.,1.),[0.,1.])
        with self.assertRaises(ValueError):eligibility(2.,1.,1.)

    def test_reused_mean_loss_matches_declared_eq7_convention(self):
        for multiple in (1,2):
            m=OrbitAveragedModel(225.,774.,.1,multiple*774e-9,0.,0.,258.,'source Eq7',charge_magnitude=0.,include_photons=False)
            r=m.rates(40.)
            self.assertAlmostEqual(r['glueball_power_gev_per_ps']*r['full_signed_period_ps']/(multiple*.1*774e-9),1.,places=13)

    def test_neutral_cooling_scaling_and_zero_emission(self):
        def model(epsilon):return OrbitAveragedModel(225.,774.,epsilon,774e-9,0.,0.,258.,'source Eq7',charge_magnitude=0.,include_photons=False)
        self.assertTrue(math.isinf(linear_cooling_time_ps(model(0.),40.,8*774e-9)))
        a=linear_cooling_time_ps(model(.1),40.,8*774e-9)
        b=linear_cooling_time_ps(model(.01),40.,8*774e-9)
        self.assertAlmostEqual(b/a,10.,places=11)


if __name__=='__main__':unittest.main()
