"""Massive endpoint flux in temporal conformal gauge (flat spacetime).

The incoming unit characteristic points into the retained string at each end:
A=a'(t+s) at the left, B=b'(t-s) at the right, for
x/c = (a(t+s)+b(t-s))/2. See studies/CONFORMAL_BOUNDARY_DERIVATION.md.
This module contains no central-force, mesh, or radial-reference imports.
"""
import numpy as np


def massive_boundary(momentum, mass, incoming):
    """Return inward label speed q, outgoing unit ray and dp/dt / tension.

    ds_left/dt=q_left, ds_right/dt=-q_right. Units: p,m in GeV;
    t,s in ps. Multiply force_per_tension by Lambda**2/hbar in GeV/ps.
    Endpoint energy transfer / tension equals q. No massless limit supplied.
    """
    p = np.asarray(momentum, dtype=float)
    k = np.asarray(incoming, dtype=float)
    if (p.shape != (3,) or k.shape != (3,) or not np.all(np.isfinite(p))
            or not np.all(np.isfinite(k)) or not np.isfinite(mass) or mass <= 0):
        raise ValueError('finite vectors and strictly positive endpoint mass required')
    if abs(np.linalg.norm(k)-1) > 1e-10:
        raise ValueError('incoming characteristic must be unit length')
    energy = np.hypot(mass, np.linalg.norm(p))
    v = p/energy
    vk = float(v@k)
    q = (vk-float(v@v))/(1-vk)
    if not -1 < q < 1:
        raise ValueError('endpoint label became null; no regular massive boundary')
    outgoing = (2*v-(1+q)*k)/(1-q)
    force = (1+q)*k-v
    return q, outgoing, force
