import math
import unittest
import numpy as np
from quirksim.pwave_annihilation import derivative_kernel,j_averages,p_projectors
from quirksim.quirk_g_parity import invariant_checks,neutral_capture_channels


class PWaveTest(unittest.TestCase):
    def test_two_color_twisted_charge_conjugation_and_channels(self):
        r=invariant_checks()
        for key in ('gauge_conjugation_residual','double_conjugation_residual','mass_splitting_residual'):self.assertLess(r[key],1e-14)
        c=neutral_capture_channels();self.assertFalse(c['3S1']['allowed_from_source_with_G_even_radiation'])
        self.assertTrue(all(c[s]['allowed_from_source_with_G_even_radiation'] for s in ('1S0','3P0','3P1','3P2')))

    def test_massless_Born_P_wave_widths_and_Landau_Yang(self):
        k=derivative_kernel(1.,0.,.01,n_theta=8,n_phi=8);j=j_averages(k)
        self.assertAlmostEqual(j['0']*3/(4*math.pi)/.01**2,9.,places=5)
        self.assertAlmostEqual(j['2']*3/(4*math.pi)/.01**2,12/5,places=5)
        self.assertLess(abs(j['1']),1e-15)
        self.assertAlmostEqual(j['2']/j['0'],4/15,places=7)
        self.assertLess(np.max(abs(sum(p_projectors().values())-np.eye(9))),1e-14)

    def test_P_projection_equal_ordered_diagrams(self):
        kwargs=dict(mass=225.,vector_mass=80.377,alpha=.03,identical=False,n_theta=8,n_phi=8)
        full=derivative_kernel(**kwargs);t=derivative_kernel(**kwargs,diagram='t');u=derivative_kernel(**kwargs,diagram='u')
        np.testing.assert_allclose(t,full/4,rtol=1e-6,atol=1e-18)
        np.testing.assert_allclose(u,full/4,rtol=1e-6,atol=1e-18)
