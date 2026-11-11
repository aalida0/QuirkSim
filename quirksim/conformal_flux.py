"""Experimental implicit boundary update using common characteristic fluxes.

For each endpoint: dP = tension/2*(incoming-outgoing) + magnetic impulse,
dX/c = (incoming+outgoing)/2. Energy transfer sets the moving label.
Outgoing histories use the same spherical segment primitive as storage.
This is a discrete balance law, not a projection of endpoint positions.
"""
import numpy as np
from scipy.optimize import root
from .conformal_boundary import massive_boundary
from .dynamics import C_CM_PER_PS as C
from .field_dynamics import MAGNETIC_FORCE_GEV_PER_PS_T


class FluxSolveFailure(RuntimeError):
    pass


def flux_trial(solver, t, y, h, caps):
    new = y.copy()
    field_impulses, displacements = [], []
    for i, incoming, outgoing in ((0, solver.a, solver.b), (1, solver.b, solver.a)):
        sign = 1. if i == 0 else -1.
        p0 = y[8+3*i:11+3*i]
        mass = solver.masses[i]
        e0 = np.hypot(mass, np.linalg.norm(p0))
        z0 = t+sign*y[i]
        cap = caps[i]
        scale = solver.force*h
        def evaluate(u):
            dp = scale*u[:3]
            p = p0+dp
            e = np.hypot(mass, np.linalg.norm(p))
            # Stable difference of energies, avoiding subtraction of ~100 GeV.
            de = (2*p0@dp+dp@dp)/(e+e0)
            # Parameterize the finite label displacement inside the causal cone.
            # Newton probes must not abort merely for trying an unphysical label.
            q = np.tanh(u[3])
            ds = h*q
            zend = z0+h+ds
            near = cap is not None and zend >= cap-32*np.spacing(max(1., abs(cap)))
            incident = incoming.left_value(cap) if near else incoming.value(zend)
            _, out, _ = massive_boundary(p, mass, incident)
            inc = incoming.primitive(zend)-incoming.primitive(z0)
            out_integral = (h-ds)*outgoing._arc(outgoing.n[-1], out, 1.)[1]
            dx_over_c = .5*(inc+out_integral)
            magnetic = (solver.charges[i]*MAGNETIC_FORCE_GEV_PER_PS_T
                        *np.cross(dx_over_c, solver.magnetic_field))
            residual = dp-.5*solver.force*(inc-out_integral)-magnetic
            return np.r_[residual/scale, de/scale-q], dp, ds, dx_over_c, magnetic, e
        incident0 = incoming.value(z0)
        q0, _, f = massive_boundary(p0, mass, incident0)
        magnetic0 = solver.charges[i]*MAGNETIC_FORCE_GEV_PER_PS_T*np.cross(p0/e0, solver.magnetic_field)
        guess = np.r_[f+magnetic0/solver.force, np.arctanh(q0)]
        try:
            solution = root(lambda u: evaluate(u)[0], guess, options={'xtol': 1e-9})
            residual, dp, ds, dx_over_c, magnetic, e = evaluate(solution.x)
        except ValueError as error:
            raise FluxSolveFailure(str(error)) from error
        # Check absolute momentum balance. Failed solves
        # stop explicitly; no fallback to an unconverged trajectory is accepted.
        if np.max(np.abs(residual))*scale > 2e-12:
            raise FluxSolveFailure('implicit flux residual exceeds 2e-12 GeV')
        new[i] = y[i]+sign*ds
        new[8+3*i:11+3*i] = p0+dp
        new[2+3*i:5+3*i] = y[2+3*i:5+3*i]+C*dx_over_c
        # Three-point Gauss integration of the interpolated momentum path.
        # This reproduces constant-force radial proper time to quadrature accuracy;
        # the general boundary scheme remains second order, not sixth order.
        nodes = (.5-np.sqrt(3/5)/2, .5, .5+np.sqrt(3/5)/2)
        weights = (5/18, 4/9, 5/18)
        new[14+i] = y[14+i]+h*mass*sum(
            weight/np.hypot(mass, np.linalg.norm(p0+node*dp))
            for node, weight in zip(nodes, weights))
        field_impulses.append(magnetic)
        displacements.append(C*dx_over_c)
    impulse = sum(field_impulses)
    new[16:19] += impulse
    for i in (0, 1):
        midpoint = y[2+3*i:5+3*i]+.5*displacements[i]
        new[19:22] += np.cross(midpoint, field_impulses[i])
    new[22:25] -= C*(t+h/2)*impulse
    if solver.track_photons:
        new[25:27] += .5*h*(solver._rhs(t, y, caps)[25:27]
                                   +solver._rhs(t+h, new, caps)[25:27])
    return new
