"""Boundary signs and fluxes, independent of trajectory integrators."""
import numpy as np
import pytest
from quirksim.conformal_boundary import massive_boundary


def test_random_massive_flux_preserves_null_ray_and_energy_balance():
    rng = np.random.default_rng(813)
    for _ in range(200):
        k = rng.normal(size=3); k /= np.linalg.norm(k)
        p = 100*rng.normal(size=3)
        v = p/np.hypot(225., np.linalg.norm(p))
        q, out, force = massive_boundary(p, 225., k)
        assert np.linalg.norm(out) == pytest.approx(1., abs=3e-15)
        np.testing.assert_allclose(((1+q)*k+(1-q)*out)/2, v, atol=3e-15)
        assert v@force == pytest.approx(q, abs=1e-15)


@pytest.mark.parametrize('beta', [-.8, -.2, 0., .4, .9])
def test_radial_force_and_label_motion(beta):
    p = 225*beta/np.sqrt(1-beta**2)
    q, out, force = massive_boundary([p, 0, 0], 225., [1., 0, 0])
    assert q == pytest.approx(beta)
    np.testing.assert_allclose(out, [-1., 0, 0], atol=2e-14)
    np.testing.assert_allclose(force, [1., 0, 0], atol=2e-14)


@pytest.mark.parametrize('mass,k', [(0., [1, 0, 0]), (225., [2, 0, 0]),
                                   (float('nan'), [1, 0, 0])])
def test_invalid_boundary_is_rejected(mass, k):
    with pytest.raises(ValueError):
        massive_boundary([0, 0, 0], mass, k)
