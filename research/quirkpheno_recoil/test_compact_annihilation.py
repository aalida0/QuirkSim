import math
import unittest
from quirksim.radiation_annihilation import (
    ALPHA_ZERO, coulomb_ground_reference, spin_singlet_two_boson_widths)


class CompactAnnihilationTest(unittest.TestCase):
    def test_parapositronium_abelian_normalization(self):
        # Remove IC interaction and divide its spectator singlet color factor.
        # Coulomb positronium gives Gamma_pPs = alpha^5 m_e / 2 at LO.
        m=.00051099895
        psi=coulomb_ground_reference(m,m,ALPHA_ZERO)['psi_squared_gev3']
        gamma=spin_singlet_two_boson_widths(m,psi,2,0.)['partial_widths_gev']['gamma_gamma']/2
        self.assertAlmostEqual(gamma/(ALPHA_ZERO**5*m/2),1.,places=13)

    def test_eta_quarkonium_color_and_charge_ratio(self):
        r=spin_singlet_two_boson_widths(10.,1.,3,.1,charge=2/3)['partial_widths_gev']
        # Standard eta_Q LO ratio for Q=2/3: 9 alpha_s^2 / (8 alpha^2).
        self.assertAlmostEqual((r['ic_gluon_pair']/r['gamma_gamma'])/(9*.1**2/(8*ALPHA_ZERO**2)),1.,places=13)
        self.assertAlmostEqual(r['ic_gluon_pair']/(8*math.pi*.1**2/(3*10**2)),1.,places=13)

    def test_unknown_total_and_formation_are_preserved(self):
        r=spin_singlet_two_boson_widths(225.,12.,2,.04)
        self.assertIsNone(r['total_width_gev']); self.assertIsNone(r['formation_probability'])
        self.assertIsNone(r['branching_fractions'])
        with self.assertRaises(ValueError):spin_singlet_two_boson_widths(225.,-1.,2,.04)
        with self.assertRaises(ValueError):spin_singlet_two_boson_widths(225.,1.,2,float('nan'))


if __name__=='__main__':unittest.main()
