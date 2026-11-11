"""Numerical radial boundary-flux check with EXPLICIT signed continuation.

This is a symmetry reduction, not a general characteristic chart continuation.
At zero length reverse the orientation of the ideal radial string, preserving
endpoint identities and momenta. This assumption supplies new chart data.
No annihilation or physical survival probability follows from this convention.
"""
import numpy as np
from scipy.integrate import solve_ivp
from .conformal_boundary import massive_boundary
from .physics import HBAR_GEV_PS
from .dynamics import C_CM_PER_PS as C


def signed_radial_boundary(mass_gev, lambda_ev, initial_length_cm,
                           initial_beta, duration_ps, max_step_ps=.5,
                           rtol=1e-10, atol=1e-12):
    """Numerically integrate r_R,p_R,tau_R; log every prescribed orientation flip."""
    if (not all(np.isfinite([mass_gev, lambda_ev, initial_length_cm,
                            initial_beta, duration_ps, max_step_ps, rtol, atol]))
            or min(mass_gev, lambda_ev, initial_length_cm, duration_ps,
                   max_step_ps, rtol, atol) <= 0 or not 0 <= initial_beta < 1):
        raise ValueError('finite positive scales and subluminal initial beta required')
    force = (lambda_ev*1e-9)**2/HBAR_GEV_PS
    momentum = mass_gev*initial_beta/np.sqrt(1-initial_beta**2)
    y = np.array([initial_length_cm/2, momentum, 0.])
    orientation, time = 1, 0.
    segments, crossings = [], []
    while time < duration_ps:
        def rhs(t, z):
            energy = np.hypot(mass_gev, z[1])
            _, _, flux = massive_boundary([z[1], 0, 0], mass_gev,
                                          [-float(orientation), 0, 0])
            return [C*z[1]/energy, force*flux[0], mass_gev/energy]

        def contact(t, z):
            return z[0]

        contact.terminal = True
        contact.direction = -orientation
        sol = solve_ivp(rhs, (time, duration_ps), y, method='DOP853',
                        max_step=max_step_ps, rtol=rtol, atol=atol,
                        dense_output=True, events=contact)
        if not sol.success:
            raise RuntimeError(sol.message)
        segments.append(sol)
        new_time = float(sol.t[-1])
        if new_time <= time:
            raise RuntimeError('radial contact restart made no progress')
        time, y = new_time, sol.y[:, -1].copy()
        if len(sol.t_events[0]):
            crossings.append({'time_ps': time, 'momentum_right_gev': float(y[1]),
                              'orientation_before': orientation,
                              'orientation_after': -orientation,
                              'prescription': 'SIGNED_RADIAL_CONTINUATION_ONLY'})
            orientation *= -1
    return segments, crossings
