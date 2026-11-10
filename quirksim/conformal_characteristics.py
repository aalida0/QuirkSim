"""Research characteristic NG integrator with moving massive boundaries.

Bulk rays propagate exactly in their null arguments; stored ray directions
use spherical interpolation, independently of the finite-element solver.
Fixed-step RK4 advances the boundary ODE, restricted by causal history.
Stops before the energy-coordinate interval collapses. No crossing/reconnection,
weak decay, radiation, external field, or physical survival prescription.
"""
from bisect import bisect_right
import numpy as np
from .conformal_boundary import massive_boundary
from .physics import HBAR_GEV_PS
from .dynamics import C_CM_PER_PS as C


class UnitRayHistory:
    """Causal unit-vector samples with exact great-circle segment primitives."""

    def __init__(self, coordinates, directions):
        self.z = []
        self.n = []
        self.integrals = []
        for z, n in zip(coordinates, directions, strict=True):
            self.append(z, n)
        if len(self.z) < 2:
            raise ValueError('at least two characteristic samples required')

    @staticmethod
    def _arc(a, b, r):
        dot = np.clip(a@b, -1., 1.)
        tangent = b-dot*a
        sine = np.linalg.norm(tangent)
        theta = np.arctan2(sine, dot)
        if theta > np.pi-1e-6:
            raise ValueError('unresolved antipodal ray; refine characteristic sampling')
        if theta < 1e-12:
            return a+r*(b-a), r*a+.5*r*r*(b-a)
        tangent /= sine
        angle = r*theta
        direction = np.cos(angle)*a+np.sin(angle)*tangent
        primitive = (np.sin(angle)*a+2*np.sin(angle/2)**2*tangent)/theta
        return direction, primitive

    def append(self, z, n):
        z = float(z)
        n = np.array(n, dtype=float, copy=True)
        if (not np.isfinite(z) or n.shape != (3,) or not np.all(np.isfinite(n))
                or abs(np.linalg.norm(n)-1) > 1e-9):
            raise ValueError('finite coordinate and unit characteristic required')
        # Remove roundoff in an analytically unit outgoing ray. Reject a physical
        # constraint violation above; repeated reflections must not amplify it.
        n /= np.linalg.norm(n)
        if self.z and z <= self.z[-1]:
            raise ValueError('characteristic coordinates must strictly increase')
        integral = np.zeros(3) if not self.z else self.integrals[-1]+(
            z-self.z[-1])*self._arc(self.n[-1], n, 1.)[1]
        self.z.append(z)
        self.n.append(n)
        self.integrals.append(integral)

    def _locate(self, z):
        tol = 1e-12*max(1., abs(z), abs(self.z[-1]))
        if z < self.z[0]-tol or z > self.z[-1]+tol:
            raise ValueError('unavailable causal characteristic history; reduce step')
        z = float(np.clip(z, self.z[0], self.z[-1]))
        i = min(max(0, bisect_right(self.z, z)-1), len(self.z)-2)
        h = self.z[i+1]-self.z[i]
        return i, h, (z-self.z[i])/h

    def value(self, z):
        i, _, r = self._locate(z)
        return self._arc(self.n[i], self.n[i+1], r)[0]

    def primitive(self, z):
        i, h, r = self._locate(z)
        return self.integrals[i]+h*self._arc(self.n[i], self.n[i+1], r)[1]


class ConformalString:
    """Two massive boundaries coupled by causal characteristic histories.

    Initial s grid is an energy coordinate in ps, not spatial arclength.
    A,B arrays give unit characteristics at t=0 on the same s grid.
    Endpoint ordering is by worldsheet label; no geometric sorting is done.
    """

    def __init__(self, mass_gev, lambda_ev, s_grid_ps, a_rays, b_rays,
                 endpoint_momenta_gev, left_position_cm=(0., 0., 0.)):
        self.masses = np.asarray(mass_gev, float)
        s = np.asarray(s_grid_ps, float)
        a = np.asarray(a_rays, float)
        b = np.asarray(b_rays, float)
        p = np.asarray(endpoint_momenta_gev, float)
        if (self.masses.shape != (2,) or not np.all(np.isfinite(self.masses))
                or np.any(self.masses <= 0) or not np.isfinite(lambda_ev)
                or lambda_ev <= 0 or s.ndim != 1 or len(s) < 2
                or not np.all(np.isfinite(s)) or np.any(np.diff(s) <= 0)
                or a.shape != (len(s), 3) or b.shape != a.shape
                or p.shape != (2, 3) or not np.all(np.isfinite(p))):
            raise ValueError('invalid conformal initial data')
        self.force = (lambda_ev*1e-9)**2/HBAR_GEV_PS
        self.a = UnitRayHistory(s, a)
        self.b = UnitRayHistory(-s[::-1], b[::-1])
        # Compatibility: endpoint worldline must lie on the initial worldsheet.
        for i, incoming, outgoing in [(0, a[0], b[0]), (1, b[-1], a[-1])]:
            _, computed, _ = massive_boundary(p[i], self.masses[i], incoming)
            if np.linalg.norm(computed-outgoing) > 1e-9:
                raise ValueError('endpoint velocity incompatible with initial rays')
        left = np.asarray(left_position_cm, float)
        if left.shape != (3,) or not np.all(np.isfinite(left)):
            raise ValueError('invalid initial position')
        chord = C/2*(self.a.primitive(s[-1])-self.a.primitive(s[0])
                     -self.b.primitive(-s[0])+self.b.primitive(-s[-1]))
        # y: labels, endpoint positions, momenta, endpoint proper times.
        self.y = np.r_[s[0], s[-1], left, left+chord, p.ravel(), 0., 0.]
        self.t = 0.

    @classmethod
    def radial(cls, mass_gev, lambda_ev, length_cm=.1, beta=.4):
        if (not np.isfinite(length_cm) or length_cm <= 0
                or not np.isfinite(beta) or not 0 <= beta < 1):
            raise ValueError('positive length and subluminal outward beta required')
        s = np.array([0., length_cm/C])
        a = np.tile([1., 0., 0.], (2, 1))
        p = mass_gev*beta/np.sqrt(1-beta**2)
        return cls([mass_gev]*2, lambda_ev, s, a, -a,
                   [[-p, 0., 0.], [p, 0., 0.]], [-length_cm/2, 0., 0.])

    def _rhs(self, t, y):
        p = y[8:14].reshape(2, 3)
        e = np.hypot(self.masses, np.linalg.norm(p, axis=1))
        v = p/e[:, None]
        ka = self.a.value(t+y[0])
        kb = self.b.value(t-y[1])
        ql, _, fl = massive_boundary(p[0], self.masses[0], ka)
        qr, _, fr = massive_boundary(p[1], self.masses[1], kb)
        return np.r_[ql, -qr, C*v.ravel(), self.force*np.r_[fl, fr],
                     self.masses/e]

    def advance(self, target_ps, max_step_ps=.02, minimum_interval_ps=1e-8):
        """Advance until target or explicit chart-collapse cutoff; never cross it."""
        if (not np.isfinite(target_ps) or target_ps < self.t
                or not np.isfinite(max_step_ps) or max_step_ps <= 0
                or not np.isfinite(minimum_interval_ps) or minimum_interval_ps <= 0):
            raise ValueError('invalid time/step/cutoff')
        while self.t < target_ps:
            interval = self.y[1]-self.y[0]
            if interval <= minimum_interval_ps:
                return 'STOP_ENERGY_COORDINATE_COLLAPSE'
            # The interval cannot shrink by more than 2h. Avoid crossing cutoff.
            h = min(max_step_ps, target_ps-self.t, .2*interval)
            if self.t+h == self.t:
                return 'STOP_FLOATING_POINT_STEP'
            t, y = self.t, self.y
            k1 = self._rhs(t, y)
            k2 = self._rhs(t+h/2, y+h*k1/2)
            k3 = self._rhs(t+h/2, y+h*k2/2)
            k4 = self._rhs(t+h, y+h*k3)
            new = y+h*(k1+2*k2+2*k3+k4)/6
            tn = t+h
            ql, bout, _ = massive_boundary(new[8:11], self.masses[0],
                                           self.a.value(tn+new[0]))
            qr, aout, _ = massive_boundary(new[11:14], self.masses[1],
                                           self.b.value(tn-new[1]))
            self.a.append(tn+new[1], aout)
            self.b.append(tn-new[0], bout)
            self.t, self.y = tn, new
        return 'TARGET_REACHED'

    def snapshot(self, samples=33):
        """Independent history integrals give string charges and chord closure."""
        if samples < 2:
            raise ValueError('at least two samples required')
        t, y = self.t, self.y
        sl, sr = y[:2]
        ia = self.a.primitive(t+sr)-self.a.primitive(t+sl)
        ib = self.b.primitive(t-sl)-self.b.primitive(t-sr)
        p = y[8:14].reshape(2, 3)
        chord = C/2*(ia-ib)
        string_momentum = self.force/2*(ia+ib)
        norm_error = 0.
        orth_error = 0.
        for s in np.linspace(sl, sr, samples):
            a, b = self.a.value(t+s), self.b.value(t-s)
            v, tangent = (a+b)/2, (a-b)/2
            norm_error = max(norm_error, abs(v@v+tangent@tangent-1))
            orth_error = max(orth_error, abs(v@tangent))
        return {
            'time_ps': t,
            'positions_cm': y[2:8].reshape(2, 3).tolist(),
            'momenta_gev': p.tolist(),
            'proper_times_ps': y[14:16].tolist(),
            'string_energy_gev': self.force*(sr-sl),
            'total_energy_gev': float(np.hypot(self.masses,
                np.linalg.norm(p, axis=1)).sum()+self.force*(sr-sl)),
            'total_momentum_gev': (p.sum(axis=0)+string_momentum).tolist(),
            'chord_closure_error_cm': float(np.linalg.norm(y[5:8]-y[2:5]-chord)),
            'max_conformal_norm_residual': float(norm_error),
            'max_conformal_orthogonality_residual': float(orth_error),
            'energy_coordinate_interval_ps': float(sr-sl),
        }
