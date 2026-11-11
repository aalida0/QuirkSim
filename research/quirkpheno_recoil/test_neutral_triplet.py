import math
import unittest
from unittest.mock import patch
import numpy as np
from quirksim.neutral_triplet import hard_spin_zero, potential, ew_inputs, ground_state, widths_from_state, gaussian_overlap_log_ratio


def two_vector_spin_kernel(m,m1,m2,alpha,identical=False):
    # Independent threshold Dirac spinors and explicit physical polarizations.
    ident=np.eye(2);zero=np.zeros((2,2))
    pa=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    gamma=[np.block([[ident,zero],[zero,-ident]])]+[np.block([[zero,s],[-s,zero]]) for s in pa]
    def slash(p):return gamma[0]*p[0]-sum((gamma[j+1]*p[j+1] for j in range(3)),np.zeros((4,4),complex))
    e1=(4*m*m+m1*m1-m2*m2)/(4*m);e2=2*m-e1;k=math.sqrt(e1*e1-m1*m1)
    p=np.array([m,0,0,0]);k1=np.array([e1,0,0,k]);k2=np.array([e2,0,0,-k])
    pol1=[np.array([0,1,0,0]),np.array([0,0,1,0])];pol2=pol1.copy()
    if m1:pol1+=[np.array([k/m1,0,0,e1/m1])]
    if m2:pol2+=[np.array([k/m2,0,0,-e2/m2])]
    denominator=(m1*m1+m2*m2)/2-2*m*m
    gram=np.zeros((4,4),complex)
    for eps1 in pol1:
        for eps2 in pol2:
            op=(slash(eps2)@(slash(p-k1)+m*np.eye(4))@slash(eps1)+slash(eps1)@(slash(p-k2)+m*np.eye(4))@slash(eps2))/denominator*(4*math.pi*alpha)
            amp=[]
            for a in range(2):
                for b in range(2):
                    u=np.eye(4)[:,b]*math.sqrt(2*m);v=np.eye(4)[:,2+a]*math.sqrt(2*m)
                    amp.append(v.conj()@gamma[0]@op@u)
            gram+=np.outer(np.conjugate(amp),amp)
    return gram*(k/m)/(32*math.pi*m*m)/(2 if identical else 1)


class NeutralTripletTest(unittest.TestCase):
    def test_radial_solver_against_analytic_coulomb_control(self):
        with patch('quirksim.neutral_triplet.ew_inputs',return_value=(0.,.25,0.)):
            state=ground_state(225.,1.,2,.04,0.,points=4000)
        self.assertLess(abs(state['energy_relative_to_neutral_threshold_gev']/(-225*.03**2/4)-1),1e-4)
        exact=(225*.03/2)**3/math.pi
        self.assertLess(abs(state['psi_at_origin_gev_three_halves'][1]**2/exact-1),1e-4)

    def test_explicit_dirac_spin_sums_and_massive_polarizations(self):
        from quirksim.physics import MW,MZ
        for m in (100.,225.,1000.):
            aw,s,ae=ew_inputs();n=2;k=hard_spin_zero(m,n,.04)
            for label,m1,m2,alpha,ident in [
                ('W_plus_W_minus',MW,MW,aw,False),('gamma_gamma',0,0,ae,True),
                ('gamma_Z',0,MZ,math.sqrt(ae*aw*(1-s)),False),('Z_Z',MZ,MZ,aw*(1-s),True)]:
                gram=two_vector_spin_kernel(m,m1,m2,alpha,ident)
                ev=np.linalg.eigvalsh(gram)
                self.assertLess(max(abs(ev[:-1]))/ev[-1],1e-12)
                index=1 if label=='W_plus_W_minus' else 0
                self.assertAlmostEqual(n*ev[-1]/k[label][index,index],1.,places=11)

    def test_coherent_flavor_interference(self):
        kernels=hard_spin_zero(225.,2,.04)
        for k in kernels.values():self.assertGreaterEqual(np.linalg.eigvalsh(k).min(),-1e-20)
        a=widths_from_state([1,0,-1],kernels)
        self.assertTrue(all(abs(v)<1e-20 for v in a.values()))
        b=widths_from_state([1,0,1],kernels)
        self.assertGreater(b['gamma_gamma'],0)

    def test_triplet_generator_potential_short_distance_spectrum(self):
        aw,s,ae=ew_inputs();r=1e-10
        v=potential([r],0.,2,.04,0.)[0]*r
        np.testing.assert_allclose(np.linalg.eigvalsh(v),-.03+aw*np.array([-2,-1,1]),atol=1e-9)

    def test_ground_state_normalization_and_flavor_symmetry(self):
        a=ground_state(225.,.157170276,2,.04,774.,points=2000)
        self.assertAlmostEqual(sum(a['flavor_probabilities']),1.,places=12)
        self.assertAlmostEqual(a['psi_at_origin_gev_three_halves'][0],a['psi_at_origin_gev_three_halves'][2],places=10)
        self.assertIsNone(a['occupation_in_transport'])

    def test_classical_means_do_not_determine_packet_overlap(self):
        self.assertLess(gaussian_overlap_log_ratio(4.1544,.001),-1e6)
        self.assertGreater(gaussian_overlap_log_ratio(4.1544,10.),-.1)


if __name__=='__main__':unittest.main()
