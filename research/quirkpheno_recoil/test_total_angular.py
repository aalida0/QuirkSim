import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.audit_triplet_total_angular import rest_angular,boost_z,HBARC,QB


class TotalAngularTest(unittest.TestCase):
    def test_point_particle_has_zero_intrinsic_orbital_action(self):
        p=np.array([5.,1.,2.,1.]);x=np.array([3.,2.,-1.,4.])
        tensor=np.outer(x,p)-np.outer(p,x)
        np.testing.assert_allclose(rest_angular(tensor,p)*HBARC,0.,atol=1e-14)

    def test_longitudinal_boost_preserves_rest_rotation(self):
        tensor=np.zeros((4,4));tensor[1,2]=7*HBARC;tensor[2,1]=-7*HBARC
        transform=np.column_stack([boost_z(np.eye(4)[i],-.4) for i in range(4)])
        p=transform@np.array([4.,0.,0.,0.])
        np.testing.assert_allclose(rest_angular(transform@tensor@transform.T,p),[0,0,7],atol=1e-12)

    def test_uniform_B_axial_impulse_matches_exact_helix(self):
        # Analytic Lorentz orbit from the origin: x=R(1-cos a), y=R sin a.
        b=2.;r=3.;a=.7;k=QB*b
        x=np.array([r*(1-np.cos(a)),r*np.sin(a),0.])
        p=np.array([k*r*np.sin(a),k*r*np.cos(a),0.])
        self.assertAlmostEqual(np.cross(x,p)[2],-k*(x@x)/2,places=15)
