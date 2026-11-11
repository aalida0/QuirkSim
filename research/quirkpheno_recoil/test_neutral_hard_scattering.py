import math
import unittest
import numpy as np
from quirksim.neutral_hard_scattering import equal_vector_kernel,bare_neutral_kernels,spin_averages
from quirksim.neutral_triplet import hard_spin_zero


class NeutralHardScatteringTest(unittest.TestCase):
    def test_qed_cross_section_from_detailed_balance(self):
        # Breit-Wheeler inverse process with spin counts and identical photons:
        # sigma_ann = sigma_BW/(2 beta^2); v_Moller=2 beta.
        for b in (.05,.4,.8):
            alpha=.01;m=225.
            k=equal_vector_kernel(m,b,0,alpha,True,n_theta=48)
            actual=spin_averages(k)['unpolarized_sigma_v_gev_minus2']
            expected=math.pi*alpha**2/(2*m*m)*(1-b*b)/b*((3-b**4)*math.log((1+b)/(1-b))-2*b*(2-b*b))
            self.assertAlmostEqual(actual/expected,1.,places=10)

    def test_neutral_threshold_projection(self):
        ref=hard_spin_zero(225.,2,.04)
        for name,k in bare_neutral_kernels(225.,0.,2,.04,n_theta=16).items():
            a=spin_averages(k)
            self.assertAlmostEqual(a['singlet_sigma_v_gev_minus2']/ref[name][1,1],1.,places=11)
            self.assertLess(abs(a['triplet_average_sigma_v_gev_minus2']),1e-19)

    def test_finite_velocity_positive_kernel_and_quadrature_convergence(self):
        a=bare_neutral_kernels(225.,.4,2,.04,n_theta=24)
        b=bare_neutral_kernels(225.,.4,2,.04,n_theta=48)
        for name,k in b.items():
            self.assertGreater(np.linalg.eigvalsh(k).min(),-1e-20)
            np.testing.assert_allclose(k,a[name],rtol=1e-9,atol=1e-18)
            self.assertGreater(spin_averages(k)['triplet_average_sigma_v_gev_minus2'],0)


if __name__=='__main__':unittest.main()
