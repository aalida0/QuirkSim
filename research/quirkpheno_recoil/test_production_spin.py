import math
import unittest
import numpy as np
from quirksim.production_spin import (GAMMA,dy_vector_density,vector_currents,singlet_fraction,
    u_spinors,v_spinors,pion_spin_map,canonical_boost_rotation,rotate_pair)
from quirksim.decays import boost,sample_decay


class ProductionSpinTest(unittest.TestCase):
    def test_vector_current_has_no_cm_singlet(self):
        for beta in (0.,.4,.9):
            p=225*beta/math.sqrt(1-beta*beta)
            for n in ([1,0,0],[0,0,1],[.6,0,.8]):
                for w in ((1,0),(0,1),(.5,.5)):
                    rho=dy_vector_density(p*np.array(n),225,helicity_weights=w)
                    self.assertLess(abs(singlet_fraction(rho)),1e-14)
                    self.assertGreater(np.linalg.eigvalsh(rho).min(),-1e-14)

    def test_vector_born_angular_distribution(self):
        beta=.4;m=225;p=m*beta/math.sqrt(1-beta*beta)
        transverse=np.sum(abs(vector_currents([p,0,0],m))**2)
        forward=np.sum(abs(vector_currents([0,0,p],m))**2)
        self.assertAlmostEqual(transverse/forward,1-beta*beta/2,places=13)

    def test_pion_eom_scalar_reduction_and_antiparticle_map(self):
        rng=np.random.default_rng(15);M=225.1571702763212;m=225.
        for p in (np.zeros(3),np.array([98.,0,63.])):
            decay=sample_decay('pion',p,M,m,rng);d=decay['daughter'];q=decay['pion']
            parent=u_spinors(p,M);daughter=u_spinors(d[1:],m)
            slash=GAMMA[0]*q[0]-sum(GAMMA[i+1]*q[i+1] for i in range(3))
            matrix=daughter.conj().T@GAMMA[0]@slash@parent
            scalar=(M-m)*(daughter.conj().T@GAMMA[0]@parent)
            np.testing.assert_allclose(matrix,scalar,rtol=1e-10,atol=1e-8)
            rotation=pion_spin_map(p,M,d[1:],m)
            anti=(v_spinors(p,M).conj().T@GAMMA[0]@v_spinors(d[1:],m)).T
            anti/=math.sqrt(np.trace(anti.conj().T@anti).real/2)
            np.testing.assert_allclose(anti,-rotation,atol=1e-13)
            if not np.any(p):np.testing.assert_allclose(rotation,np.eye(2),atol=1e-14)

    def test_wigner_roundtrip_and_frame_dependent_spin_fraction(self):
        m=225.;p=np.array([m*.4/math.sqrt(.84),0,0]);b=np.array([0,0,.25])
        rho=dy_vector_density(p,m);maps=[]
        for pp in (p,-p):
            w=canonical_boost_rotation(pp,m,b)
            lab=boost(math.hypot(m,np.linalg.norm(pp)),pp,b)
            back=canonical_boost_rotation(lab[1:],m,-b)
            np.testing.assert_allclose(back@w,np.eye(2),atol=1e-13)
            maps.append(w)
        labrho=rotate_pair(rho,*maps)
        self.assertGreater(singlet_fraction(labrho),0)
        restored=rotate_pair(labrho,maps[0].conj().T,maps[1].conj().T)
        self.assertLess(abs(singlet_fraction(restored)),1e-14)


if __name__=='__main__':unittest.main()
