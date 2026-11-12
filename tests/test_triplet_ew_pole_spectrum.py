"""Physics and numerical validation of the minimal Dirac Y=0 triplet EW splitting."""
from pathlib import Path
import math
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from triplet_ew_pole_spectrum import loop_f, triplet_spectrum

class TestEWTripletPoleMass(unittest.TestCase):
    def test_mg5_style_gmu_inputs_mass_dependent(self):
        # Independent scipy quadrature of Cirelli-Fornengo-Strumia one-loop
        # integral at alphaEW(MZ)^-1=127.9, GF=1.1663787e-5, MZ=91.1876
        val=triplet_spectrum(225.)
        self.assertAlmostEqual(val.mw_gev,79.82449037448939,places=8)
        self.assertAlmostEqual(val.delta_m_mev,161.8440131071835,places=5)
        self.assertEqual(val.neutral_pole_mass_gev,225.)
        self.assertAlmostEqual(val.charged_pole_mass_gev,225.1618440131,places=8)
        self.assertEqual(val.scheme,"GF_MZ_ALPHA_MZ_TREE")

    def test_mass_scan_and_heavy_limit(self):
        a=[triplet_spectrum(m).delta_m_mev for m in (50,100,150,225,300,1000,10000)]
        self.assertTrue(all(0 < val < 200 for val in a))
        self.assertTrue(all(a[i]<a[i+1] for i in range(len(a)-1)))
        # Alpha(MZ) input convention asymptotic formula
        obj=triplet_spectrum(10000)
        expected=obj.alpha2*obj.mw_gev*math.sin(
            math.asin(math.sqrt(obj.sin2theta_w))/2
        )**2*1000
        self.assertAlmostEqual(a[-1],expected,delta=0.07)

    def test_loop_integral(self):
        self.assertEqual(loop_f(0.),-5.)
        self.assertAlmostEqual(loop_f(1.),-0.5586019074,delta=0.5) # coarse sanity only
        self.assertTrue(math.isfinite(loop_f(0.00001)))

    def test_external_mw_changes_scheme(self):
        x=triplet_spectrum(225., mw_gev=80.379)
        y=triplet_spectrum(225.)
        self.assertEqual(x.scheme,"MW_MZ_ALPHA_INPUT")
        self.assertNotAlmostEqual(x.delta_m_mev,y.delta_m_mev,places=1)
        self.assertGreater(x.delta_m_mev,0)

    def test_reject_bad_inputs(self):
        for m in (0,-1,float("nan"),float("inf")):
            with self.assertRaises(ValueError):
                triplet_spectrum(m)
        with self.assertRaises(ValueError):
            triplet_spectrum(225.,mw_gev=92.)
        with self.assertRaises(ValueError):
            loop_f(-.3)
        with self.assertRaises(ValueError):
            triplet_spectrum(225.,aewm1=0)
        with self.assertRaises(ValueError):
            triplet_spectrum(225.,gf=1e-9)


if __name__=="__main__":
    unittest.main()
