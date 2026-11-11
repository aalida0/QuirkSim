"""Research characteristic NG integrator with moving massive boundaries.

Bulk rays propagate exactly in their null arguments; stored ray directions
use spherical interpolation, independently of the finite-element solver.
Fixed-step RK4 advances the boundary ODE, restricted by causal history.
Stops before the energy-coordinate interval collapses. Uniform magnetic fields
act on massive boundaries only; no bulk string electric charge is assumed.
"""
from bisect import bisect_right, bisect_left
import json
import time
import numpy as np
from .conformal_boundary import massive_boundary
from .physics import HBAR_GEV_PS
from .dynamics import C_CM_PER_PS as C
from .dynamics import HBARC_GEV_CM
from .field_dynamics import MAGNETIC_FORCE_GEV_PER_PS_T
from .radiation_annihilation import lienard_power_gev_per_ps


_GAUSS_NODES, _GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(5)


class UnitRayHistory:
    """Causal unit-vector samples with exact great-circle segment primitives."""

    def __init__(self, coordinates, directions):
        self.z = []
        self.n = []
        self.integrals = []
        self.jumps = []
        self.derivatives = []
        self.derivative_fronts = []
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

    @staticmethod
    def _checked_derivative(n, derivative):
        if derivative is None:
            return None
        derivative = np.asarray(derivative,float)
        if (derivative.shape != (3,) or not np.all(np.isfinite(derivative))
                or abs(n@derivative) > 1e-8*max(1.,np.linalg.norm(derivative))):
            raise ValueError('finite tangent derivative required')
        return derivative.copy()

    def append(self, z, n, derivative=None):
        z = float(z)
        n = np.array(n, dtype=float, copy=True)
        if (not np.isfinite(z) or n.shape != (3,) or not np.all(np.isfinite(n))
                or abs(np.linalg.norm(n)-1) > 1e-9):
            raise ValueError('finite coordinate and unit characteristic required')
        # Remove roundoff in an analytically unit outgoing ray. Reject a physical
        # constraint violation above; repeated reflections must not amplify it.
        n /= np.linalg.norm(n)
        derivative = self._checked_derivative(n,derivative)
        if self.z and z <= self.z[-1]:
            raise ValueError('characteristic coordinates must strictly increase')
        integral = np.zeros(3) if not self.z else self.integrals[-1]+(
            z-self.z[-1])*self._arc(self.n[-1], n, 1.)[1]
        self.z.append(z)
        self.n.append(n)
        self.integrals.append(integral)
        self.derivatives.append(None if derivative is None else np.asarray(derivative, float).copy())
        if len(self.z) > 1 and self.derivatives[-2] is not None and derivative is not None:
            self.integrals[-1] = self.integrals[-2]+self._segment_integral(len(self.z)-2, 1.)

    def _segment(self, i, r):
        """Unit direction and derivative per null coordinate on a smooth segment."""
        h = self.z[i+1]-self.z[i]
        a, b = self.n[i], self.n[i+1]
        da, db = self.derivatives[i], self.derivatives[i+1]
        if da is not None and db is not None:
            q = ((2*r**3-3*r*r+1)*a+(r**3-2*r*r+r)*h*da
                 +(-2*r**3+3*r*r)*b+(r**3-r*r)*h*db)
            dq = ((6*r*r-6*r)*a+(3*r*r-4*r+1)*h*da
                  +(-6*r*r+6*r)*b+(3*r*r-2*r)*h*db)
            norm = np.linalg.norm(q)
            if norm < .5:
                raise ValueError('unresolved cubic characteristic segment')
            n = q/norm
            return n, (dq-n*(n@dq))/(norm*h)
        n = self._arc(a,b,r)[0]
        dot = np.clip(a@b,-1.,1.)
        tangent = b-dot*a
        sine = np.linalg.norm(tangent)
        theta = np.arctan2(sine,dot)
        if theta < 1e-12:
            return n, (b-a)/h
        tangent /= sine
        return n, theta/h*(-np.sin(r*theta)*a+np.cos(r*theta)*tangent)

    def _segment_integral(self, i, r):
        h = self.z[i+1]-self.z[i]
        if self.derivatives[i] is None or self.derivatives[i+1] is None:
            return h*self._arc(self.n[i],self.n[i+1],r)[1]
        # Fixed five-point Gauss quadrature of the normalized Hermite curve.
        nodes, weights = _GAUSS_NODES, _GAUSS_WEIGHTS
        return h*r/2*sum(w*self._segment(i,r*(x+1)/2)[0] for x,w in zip(nodes,weights))

    def _locate(self, z):
        tol = 1e-12*max(1., abs(z), abs(self.z[-1]))
        if z < self.z[0]-tol or z > self.z[-1]+tol:
            raise ValueError('unavailable causal characteristic history; reduce step')
        z = float(np.clip(z, self.z[0], self.z[-1]))
        i = max(0, bisect_right(self.z, z)-1)
        if i == len(self.z)-1:
            return i, 0., 0.
        h = self.z[i+1]-self.z[i]
        return i, h, (z-self.z[i])/h

    def _jump_coordinate(self, z):
        if self.jumps:
            j = bisect_left(self.jumps, z)
            for candidate in self.jumps[max(0, j-1):j+1]:
                if abs(z-candidate) <= 16*np.spacing(max(1., abs(candidate))):
                    z = candidate
                    break
        return z

    def value(self, z):
        z = self._jump_coordinate(z)
        i, _, r = self._locate(z)
        if i == len(self.z)-1:
            return self.n[i].copy()
        if self.derivatives[i] is not None and self.derivatives[i+1] is not None:
            return self._segment(i,r)[0]
        return self._arc(self.n[i], self.n[i+1], r)[0]

    def gradient(self, z, left=False):
        if not left:
            z = self._jump_coordinate(z)
        if left:
            index = bisect_left(self.z,z)
            if index < len(self.z) and self.z[index] == z and index > 0:
                return self._segment(index-1,1.)[1]
        i, _, r = self._locate(z)
        if i == len(self.z)-1:
            if self.derivatives[i] is not None:
                return self.derivatives[i].copy()
            return self._segment(i-1,1.)[1]
        return self._segment(i,r)[1]

    def left_value(self, z):
        i = bisect_left(self.z, z)
        if i < len(self.z) and self.z[i] == z:
            return self.n[i].copy()
        return self.value(z)

    def next_jump(self, z):
        i = bisect_right(self.jumps, z+32*np.spacing(max(1., abs(z))))
        return self.jumps[i] if i < len(self.jumps) else None

    def primitive(self, z):
        i, h, r = self._locate(z)
        if h == 0:
            return self.integrals[i].copy()
        return self.integrals[i]+self._segment_integral(i,r)

    def jump(self, z, n, derivative=None, derivative_front=False):
        """One-sided outgoing ray discontinuity; its integral has no impulse.

        The left limit is retained at the same null coordinate. The right limit
        propagates causally; no interpolation across a finite interval smears it.
        """
        n = np.array(n, dtype=float, copy=True)
        if (z != self.z[-1] or n.shape != (3,) or not np.all(np.isfinite(n))
                or abs(np.linalg.norm(n)-1) > 1e-9):
            raise ValueError('jump requires the current boundary coordinate and a unit ray')
        n /= np.linalg.norm(n)
        derivative = self._checked_derivative(n,derivative)
        self.z.append(float(z))
        self.n.append(n)
        self.integrals.append(self.integrals[-1].copy())
        self.jumps.append(float(z))
        if derivative_front:
            self.derivative_fronts.append(float(z))
        elif z in self.derivative_fronts:
            # A direction jump at the same coordinate supersedes a slope front.
            self.derivative_fronts = [point for point in self.derivative_fronts if point != z]
        self.derivatives.append(None if derivative is None else np.asarray(derivative,float).copy())


class ConformalString:
    """Two massive boundaries coupled by causal characteristic histories.

    Initial s grid is an energy coordinate in ps, not spatial arclength.
    A,B arrays give unit characteristics at t=0 on the same s grid.
    Endpoint ordering is by worldsheet label; no geometric sorting is done.
    """

    def __init__(self, mass_gev, lambda_ev, s_grid_ps, a_rays, b_rays,
                 endpoint_momenta_gev, left_position_cm=(0., 0., 0.),
                 charges=(1., -1.), magnetic_field_tesla=(0., 0., 0.),
                 track_photons=False):
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
        self.charges = np.asarray(charges, float).copy()
        self.track_photons = bool(track_photons)
        self.magnetic_field = np.asarray(magnetic_field_tesla, float).copy()
        if (self.charges.shape != (2,) or self.magnetic_field.shape != (3,)
                or not np.all(np.isfinite(self.charges))
                or not np.all(np.isfinite(self.magnetic_field))):
            raise ValueError('finite endpoint charges and magnetic field required')
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
        self.y = np.r_[s[0], s[-1], left, left+chord, p.ravel(), 0., 0., np.zeros(11)]
        self.t = 0.
        self.emissions = []
        self.reflections = []
        self.last_advance_diagnostic = None
        self.cubic_history = False

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

    @classmethod
    def boosted_radial(cls, mass_gev, lambda_ev, rest_length_cm=.1,
                       beta_rel=.4, beta_cm=.25, b_tesla=2., track_photons=False):
        """Same longitudinally boosted initial slice as the FEM benchmark.

        Length is Lorentz contracted. Equal lab-time endpoint velocities are
        velocity-added; internal longitudinal labelling is gauge freedom.
        """
        if (not np.all(np.isfinite([mass_gev, lambda_ev, rest_length_cm,
                                   beta_rel, beta_cm, b_tesla]))
                or mass_gev <= 0 or lambda_ev <= 0 or rest_length_cm <= 0
                or not 0 <= beta_rel < 1 or not abs(beta_cm) < 1):
            raise ValueError('invalid boosted radial benchmark')
        length = rest_length_cm*np.sqrt(1-beta_cm**2)
        velocities = np.array([(-beta_rel+beta_cm)/(1-beta_rel*beta_cm),
                               (beta_rel+beta_cm)/(1+beta_rel*beta_cm)])
        momenta = np.zeros((2, 3))
        momenta[:, 0] = mass_gev*velocities/np.sqrt(1-velocities**2)
        a = np.tile([1., 0., 0.], (2, 1))
        return cls([mass_gev]*2, lambda_ev, [0., length/C], a, -a,
                   momenta, [-length/2, 0., 0.], magnetic_field_tesla=[0., 0., b_tesla],
                   track_photons=track_photons)

    def boundary_ray_derivative(self, endpoint, t, y, cap=None, mass=None, charge=None, right_cap=None):
        history = self.a if endpoint == 0 else self.b
        z = t+(y[0] if endpoint == 0 else -y[1])
        left = cap is not None and z >= cap-32*np.spacing(max(1.,abs(cap)))
        if left: z = cap
        if right_cap is not None:
            z, left = right_cap, False
        k = history.left_value(z) if left else history.value(z)
        kz = history.gradient(z,left=left)
        p = y[8+3*endpoint:11+3*endpoint]
        m = self.masses[endpoint] if mass is None else mass
        charge = self.charges[endpoint] if charge is None else charge
        energy = np.hypot(m,np.linalg.norm(p))
        v = p/energy
        q, outgoing, force = massive_boundary(p,m,k)
        force = self.force*force+charge*MAGNETIC_FORCE_GEV_PER_PS_T*np.cross(v,self.magnetic_field)
        dv = (force-v*(v@force))/energy
        dk = (1+q)*kz
        vk = v@k; vv=v@v
        d_vk = dv@k+v@dk
        dq = ((d_vk-2*v@dv)*(1-vk)+(vk-vv)*d_vk)/(1-vk)**2
        slope = (2*dv-dq*k-(1+q)*dk+dq*outgoing)/(1-q)**2
        return slope-outgoing*(outgoing@slope)

    def enable_cubic_history(self):
        if self.cubic_history:
            return
        self.cubic_history = True
        da = self.boundary_ray_derivative(1,self.t,self.y)
        db = self.boundary_ray_derivative(0,self.t,self.y)
        # Initial bulk rays and newly emitted rays can agree in direction while
        # having different slopes. Carry that front causally, like a weak kink.
        self.a.jump(self.a.z[-1],self.a.n[-1],da,derivative_front=True)
        self.b.jump(self.b.z[-1],self.b.n[-1],db,derivative_front=True)

    def _rhs(self, t, y, incoming_caps=(None, None)):
        p = y[8:14].reshape(2, 3)
        e = np.hypot(self.masses, np.linalg.norm(p, axis=1))
        v = p/e[:, None]
        def incoming(history, z, cap):
            near = cap is not None and z >= cap-32*np.spacing(max(1., abs(cap)))
            return history.left_value(cap) if near else history.value(z)
        ka = incoming(self.a, t+y[0], incoming_caps[0])
        kb = incoming(self.b, t-y[1], incoming_caps[1])
        ql, _, fl = massive_boundary(p[0], self.masses[0], ka)
        qr, _, fr = massive_boundary(p[1], self.masses[1], kb)
        magnetic = (self.charges[:, None]*MAGNETIC_FORCE_GEV_PER_PS_T
                    *np.cross(v, self.magnetic_field))
        total_force = self.force*np.array([fl, fr])+magnetic
        impulse = magnetic.sum(axis=0)
        torque = np.cross(y[2:8].reshape(2, 3), magnetic).sum(axis=0)
        powers = np.zeros(2)
        if self.track_photons:
            for i in (0, 1):
                acceleration = (total_force[i]-v[i]*(v[i]@total_force[i]))/e[i]
                powers[i] = lienard_power_gev_per_ps(v[i], acceleration, self.charges[i])
        return np.r_[ql, -qr, C*v.ravel(), total_force.ravel(),
                     self.masses/e, impulse, torque, -C*t*impulse, powers]

    def advance(self, *args, **kwargs):
        from .conformal_flux import FluxSolveFailure
        try:
            return self._advance_impl(*args, **kwargs)
        except FluxSolveFailure as error:
            self.last_advance_diagnostic = {'reason': 'STOP_FLUX_SOLVE',
                'time_ps': self.t, 'message': str(error)}
            return 'STOP_FLUX_SOLVE'

    def _advance_impl(self, target_ps, max_step_ps=.02, minimum_interval_ps=1e-8,
                causal_fraction=.05, max_ray_angle=.025, max_steps=200000,
                matching_radius_cm=None, max_rejected_ray_trials=None,
                minimum_step_ulps=0, deadline_monotonic=None, integration_scheme="rk4"):
        """Advance until target or explicit chart-collapse cutoff; never cross it."""
        if (not np.isfinite(target_ps) or target_ps < self.t
                or not np.isfinite(max_step_ps) or max_step_ps <= 0
                or not np.isfinite(minimum_interval_ps) or minimum_interval_ps <= 0
                or not 0 < causal_fraction < .25 or not 0 < max_ray_angle < .5
                or max_steps < 1 or (matching_radius_cm is not None and
                    (not np.isfinite(matching_radius_cm) or matching_radius_cm <= 0))):
            raise ValueError('invalid time/step/cutoff')
        if (max_rejected_ray_trials is not None and max_rejected_ray_trials < 1) or minimum_step_ulps < 0:
            raise ValueError('invalid diagnostic resource/resolution limit')
        if integration_scheme == 'flux' and self.cubic_history:
            raise ValueError('flux trial currently requires spherical segment history')
        if integration_scheme not in ('rk4', 'flux'):
            raise ValueError('unknown boundary integration scheme')
        steps = 0
        rejected_ray_trials = 0
        event_root_truncations = 0
        min_accepted_step_ps = None
        trial_context = "no trial yet"
        self.last_advance_diagnostic = None
        def stop(reason, state):
            # Capture only at a stop; do not modify trial steps, histories or state.
            diagnostic = {'reason': reason, 'time_ps': float(self.t),
                'target_ps': float(target_ps), 'attempted_steps': int(state['steps']),
                'max_steps': int(max_steps),
                'max_rejected_ray_trials': max_rejected_ray_trials,
                'minimum_step_ulps': minimum_step_ulps,
                'rejected_ray_trials': rejected_ray_trials,
                'event_root_truncations': event_root_truncations,
                'min_accepted_step_ps': min_accepted_step_ps,
                'null_coordinate_ulps_ps': [float(np.spacing(self.a.z[-1])),
                                            float(np.spacing(self.b.z[-1]))],
                'trial_context': trial_context,
                'energy_coordinate_interval_ps': float(self.y[1]-self.y[0]),
                'ray_history_samples': [len(self.a.z), len(self.b.z)],
                'last_outgoing_null_coordinates_ps': [self.a.z[-1], self.b.z[-1]]}
            for key in ('h', 'angle', 'angle_step_ps', 'rejections', 'az', 'bz'):
                if key in state:
                    diagnostic[key] = float(state[key])
            if 'caps' in state:
                diagnostic['incoming_jump_caps_ps'] = state['caps']
            if 'roots' in state:
                diagnostic['event_roots_ps'] = state['roots']
            if 'az' in state:
                diagnostic['proposed_null_increments_ps'] = [
                    float(state['az']-self.a.z[-1]), float(state['bz']-self.b.z[-1])]
            self.last_advance_diagnostic = diagnostic
            return reason
        while self.t < target_ps:
            steps += 1
            if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
                return stop('STOP_DIAGNOSTIC_WALL_TIME', locals())
            if steps > max_steps:
                return stop('STOP_STEP_BUDGET', locals())
            if matching_radius_cm is not None and np.linalg.norm(self.y[5:8]-self.y[2:5]) <= matching_radius_cm*(1+1e-9):
                return 'QUANTUM_MATCHING_SURFACE'
            if target_ps-self.t <= 4*np.spacing(max(1., abs(target_ps))):
                return 'TARGET_REACHED'
            interval = self.y[1]-self.y[0]
            if interval <= minimum_interval_ps:
                return stop('STOP_ENERGY_COORDINATE_COLLAPSE', locals())
            # The interval cannot shrink by more than 2h. Avoid crossing cutoff.
            h = min(max_step_ps, target_ps-self.t, causal_fraction*interval)
            if self.t+h == self.t:
                return stop('STOP_FLOATING_POINT_STEP', locals())
            t, y = self.t, self.y
            caps = [self.a.next_jump(t+y[0]), self.b.next_jump(t-y[1])]
            first_stage = self._rhs(t, y, caps)
            trial_cache = {0.: y}
            def trial(h):
                if h in trial_cache:
                    return trial_cache[h]
                if integration_scheme == 'flux':
                    from .conformal_flux import flux_trial
                    trial_cache[h] = flux_trial(self, t, y, h, caps)
                    return trial_cache[h]
                k1 = first_stage
                k2 = self._rhs(t+h/2, y+h*k1/2, caps)
                k3 = self._rhs(t+h/2, y+h*k2/2, caps)
                k4 = self._rhs(t+h, y+h*k3, caps)
                result = y+h*(k1+2*k2+2*k3+k4)/6
                trial_cache[h] = result
                return result

            def incident(history, z, cap):
                near = cap is not None and z >= cap-32*np.spacing(max(1., abs(cap)))
                return history.left_value(cap) if near else history.value(z)

            rejections = 0
            while True:
                if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
                    return stop('STOP_DIAGNOSTIC_WALL_TIME', locals())
                if minimum_step_ulps:
                    resolution = max(np.spacing(max(1., abs(t))),
                                     np.spacing(max(1., abs(self.a.z[-1]))),
                                     np.spacing(max(1., abs(self.b.z[-1]))))
                    if h <= minimum_step_ulps*resolution:
                        if h == target_ps-self.t:
                            # A sub-ULP observation remainder is not a stalled
                            # physical trial. As in the existing target test,
                            # return without changing time, coordinates or state.
                            return 'TARGET_REACHED'
                        return stop('STOP_DIAGNOSTIC_TIME_RESOLUTION', locals())
                new = trial(h)
                trial_context = "current unaccepted trial"
                # Integrate to an incoming discontinuity before reflecting it.
                # Boundary momentum is continuous; outgoing characteristic is not.
                from scipy.optimize import brentq
                roots = []
                for which, cap in enumerate(caps):
                    def argument(dt):
                        yy = trial(dt)
                        return t+dt+(yy[0] if which == 0 else -yy[1])
                    end_argument = t+h+(new[0] if which == 0 else -new[1])
                    if cap is not None and end_argument >= cap:
                        roots.append(brentq(lambda dt: argument(dt)-cap, 0., h,
                                            xtol=1e-13, rtol=1e-13))
                if roots:
                    event_root_truncations += 1
                    h = min(roots)
                    if minimum_step_ulps and h <= minimum_step_ulps*resolution:
                        return stop('STOP_DIAGNOSTIC_TIME_RESOLUTION', locals())
                    new = trial(h)
                tn = t+h
                ql, bout, _ = massive_boundary(new[8:11], self.masses[0],
                                               incident(self.a, tn+new[0], caps[0]))
                qr, aout, _ = massive_boundary(new[11:14], self.masses[1],
                                               incident(self.b, tn-new[1], caps[1]))
                angle_step_ps = h
                angle = max(np.arccos(np.clip(aout@self.a.n[-1], -1., 1.)),
                            np.arccos(np.clip(bout@self.b.n[-1], -1., 1.)))
                if angle <= max_ray_angle:
                    break
                h /= 2
                trial_context = "rejected trial; h is proposed for next trial"
                rejections += 1
                rejected_ray_trials += 1
                if max_rejected_ray_trials is not None and rejected_ray_trials >= max_rejected_ray_trials:
                    return stop('STOP_DIAGNOSTIC_REJECTION_BUDGET', locals())
                if rejections > 64:
                    return stop('STOP_UNRESOLVED_NULL_RAY', locals())
                if t+h == t:
                    return stop('STOP_FLOATING_POINT_RAY_RESOLUTION', locals())
            matched = False
            if matching_radius_cm is not None:
                from scipy.optimize import brentq, minimize_scalar
                r0 = np.linalg.norm(y[5:8]-y[2:5])
                r1 = np.linalg.norm(new[5:8]-new[2:5])
                # Endpoint speeds <=c give a rigorous no-hit lower bound.
                if min(r0, r1)-2*C*h <= matching_radius_cm:
                    def surface(dt):
                        yy = trial(dt)
                        return np.linalg.norm(yy[5:8]-yy[2:5])-matching_radius_cm
                    crossing = h
                    if r1 > matching_radius_cm:
                        closest = minimize_scalar(surface, bounds=(0., h), method='bounded',
                                                  options={'xatol': max(1e-14, h*1e-10)})
                        crossing = closest.x
                    if surface(crossing) <= 0:
                        h = brentq(surface, 0., crossing, xtol=1e-13, rtol=1e-13)
                        new, tn = trial(h), t+h
                        _, bout, _ = massive_boundary(new[8:11], self.masses[0], incident(self.a, tn+new[0], caps[0]))
                        _, aout, _ = massive_boundary(new[11:14], self.masses[1], incident(self.b, tn-new[1], caps[1]))
                        matched = True
            az, bz = tn+new[1], tn-new[0]
            if az <= self.a.z[-1] or bz <= self.b.z[-1]:
                return stop('STOP_FLOATING_POINT_CHARACTERISTIC_LABEL', locals())
            min_accepted_step_ps = h if min_accepted_step_ps is None else min(min_accepted_step_ps, h)
            da = self.boundary_ray_derivative(1,tn,new,caps[1]) if self.cubic_history else None
            db = self.boundary_ray_derivative(0,tn,new,caps[0]) if self.cubic_history else None
            self.a.append(az, aout, da)
            self.b.append(bz, bout, db)
            self.t, self.y = tn, new
            trial_context = "last accepted step"
            for which, cap in enumerate(caps):
                argument = tn+(new[0] if which == 0 else -new[1])
                if cap is not None and abs(argument-cap) <= 32*np.spacing(max(1., abs(cap))):
                    if which == 0:
                        _, reflected, _ = massive_boundary(new[8:11], self.masses[0], self.a.value(cap))
                        derivative = self.boundary_ray_derivative(0,tn,new,right_cap=cap) if self.cubic_history else None
                        front = cap in self.a.derivative_fronts
                        self.b.jump(bz, reflected, derivative,derivative_front=front)
                    else:
                        _, reflected, _ = massive_boundary(new[11:14], self.masses[1], self.b.value(cap))
                        derivative = self.boundary_ray_derivative(1,tn,new,right_cap=cap) if self.cubic_history else None
                        front = cap in self.b.derivative_fronts
                        self.a.jump(az, reflected, derivative,derivative_front=front)
                    self.reflections.append({'time_ps': tn, 'endpoint': which,
                                             'incoming_null_coordinate_ps': cap})
                    if self.cubic_history:
                        self.reflections[-1]['front_kind'] = 'derivative' if front else 'direction'
            if matched:
                return 'QUANTUM_MATCHING_SURFACE'
        return 'TARGET_REACHED'

    def save_checkpoint(self, path):
        """Save exact accumulated histories for classical replay, with no pickle."""
        payload = {key: np.asarray(getattr(self, key)) for key in
                   ('y', 'masses', 'charges', 'magnetic_field')}
        payload.update(t=np.asarray(self.t), force=np.asarray(self.force),
                       track_photons=np.asarray(self.track_photons))
        for name in ('a', 'b'):
            history = getattr(self, name)
            for field in ('z', 'n', 'integrals', 'jumps'):
                payload[name+'_'+field] = np.asarray(getattr(history, field))
        metadata = dict(format=1, emissions=self.emissions, reflections=self.reflections,
                        last_advance_diagnostic=self.last_advance_diagnostic)
        metadata['cubic_history'] = self.cubic_history
        for name in ('a','b'):
            payload[name+'_derivative_fronts'] = np.asarray(getattr(self,name).derivative_fronts)
            payload[name+'_derivatives'] = np.asarray([
                [np.nan]*3 if v is None else v for v in getattr(self,name).derivatives])
        payload['metadata'] = np.asarray(json.dumps(metadata))
        # Open explicitly so numpy never silently appends a filename suffix.
        with open(path, 'wb') as handle:
            np.savez_compressed(handle, **payload)

    @classmethod
    def load_checkpoint(cls, path):
        """Restore solver-owned checkpoints; this does not resume a weak-decay scheduler."""
        with np.load(path, allow_pickle=False) as payload:
            metadata = json.loads(str(payload['metadata']))
            if metadata['format'] != 1:
                raise ValueError('unsupported characteristic checkpoint version')
            solver = cls.__new__(cls)
            for key in ('y', 'masses', 'charges', 'magnetic_field'):
                setattr(solver, key, payload[key].copy())
            solver.t, solver.force = float(payload['t']), float(payload['force'])
            solver.track_photons = bool(payload['track_photons'])
            for name in ('a', 'b'):
                history = UnitRayHistory.__new__(UnitRayHistory)
                for field in ('z', 'jumps'):
                    setattr(history, field, payload[name+'_'+field].tolist())
                for field in ('n', 'integrals'):
                    setattr(history, field, [row.copy() for row in payload[name+'_'+field]])
                key = name+'_derivatives'
                history.derivative_fronts = payload[name+'_derivative_fronts'].tolist() if name+'_derivative_fronts' in payload else []
                history.derivatives = ([None if np.isnan(row).all() else row.copy() for row in payload[key]]
                                       if key in payload else [None]*len(history.z))
                setattr(solver, name, history)
            solver.cubic_history = metadata.get('cubic_history',False)
            solver.emissions = metadata['emissions']
            solver.reflections = metadata['reflections']
            solver.last_advance_diagnostic = metadata['last_advance_diagnostic']
        return solver

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
            'external_impulse_gev': y[16:19].tolist(),
            'external_angular_impulse_gev_cm': y[19:22].tolist(),
            'external_boost_impulse_gev_cm': y[22:25].tolist(),
            'smooth_photon_energy_without_backreaction_gev': y[25:27].tolist(),
            'reflected_discontinuities': len(self.reflections),
        }

    def poincare(self, samples=129):
        """Orbital angular/boost charges with independent composite quadrature.

        Refine samples to check quadrature; these are orbital classical string
        charges. Spin of endpoints and weak products is not tracked here.
        """
        from scipy.integrate import simpson
        if samples < 3 or samples % 2 != 1:
            raise ValueError('odd quadrature sample count >= 3 required')
        t, y = self.t, self.y
        sl, sr = y[:2]
        grid = np.linspace(sl, sr, samples)
        x, density = [], []
        aleft = self.a.primitive(t+sl)
        bleft = self.b.primitive(t-sl)
        for s in grid:
            a, b = self.a.value(t+s), self.b.value(t-s)
            position = y[2:5]+C/2*(self.a.primitive(t+s)-aleft
                                   +self.b.primitive(t-s)-bleft)
            x.append(position)
            density.append(self.force*(a+b)/2)
        x, density = np.array(x), np.array(density)
        pend = y[8:14].reshape(2, 3)
        xend = y[2:8].reshape(2, 3)
        eend = np.hypot(self.masses, np.linalg.norm(pend, axis=1))
        angular = np.cross(xend, pend).sum(axis=0)+simpson(
            np.cross(x, density), x=grid, axis=0)
        energy_moment = (eend[:, None]*xend).sum(axis=0)+self.force*simpson(x, x=grid, axis=0)
        snap = self.snapshot(samples=3)
        boost = energy_moment-C*t*np.array(snap['total_momentum_gev'])
        return {'orbital_angular_momentum_gev_cm': angular.tolist(),
                'orbital_angular_momentum_hbar': (angular/HBARC_GEV_CM).tolist(),
                'boost_charge_gev_cm': boost.tolist(),
                'orbital_angular_minus_external_gev_cm': (angular-y[19:22]).tolist(),
                'boost_minus_external_gev_cm': (boost-y[22:25]).tolist(),
                'quadrature_samples': samples}

    def apply_endpoint_decay(self, endpoint, daughter_mass_gev, products):
        """Local on-shell weak impulse with a causal outgoing-ray jump.

        products contains the daughter and emitted lab four-vectors (GeV).
        The incoming characteristic and bulk histories are unchanged. Unlike
        finite-element canonical impulses no global velocity reconstruction is
        required: this gauge's boundary particle has mechanical momentum p.
        """
        if endpoint not in (0, 1) or self.charges[endpoint] == 0:
            raise ValueError('a charged endpoint label 0 or 1 is required')
        if (not np.isfinite(daughter_mass_gev) or daughter_mass_gev <= 0
                or daughter_mass_gev >= self.masses[endpoint]):
            raise ValueError('daughter must have a smaller positive mass')
        products = {k: np.asarray(v, float) for k, v in products.items()}
        if 'daughter' not in products or len(products) < 2 or any(
                p.shape != (4,) or not np.all(np.isfinite(p)) or p[0] <= 0
                or p[0]**2-np.dot(p[1:], p[1:]) < -1e-8
                for p in products.values()):
            raise ValueError('future-directed on-shell decay four-vectors required')
        start = 8+3*endpoint
        oldp = self.y[start:start+3]
        parent = np.r_[np.hypot(self.masses[endpoint], np.linalg.norm(oldp)), oldp]
        parent_charge = float(self.charges[endpoint])
        daughter = products['daughter']
        if abs(daughter[0]-np.hypot(daughter_mass_gev, np.linalg.norm(daughter[1:]))) > 1e-9:
            raise ValueError('daughter mass-shell mismatch')
        residual = parent-sum(products.values(), np.zeros(4))
        if np.max(np.abs(residual)) > 1e-8:
            raise ValueError('weak impulse violates four-momentum conservation')
        incoming = (self.a.value(self.t+self.y[0]) if endpoint == 0
                    else self.b.value(self.t-self.y[1]))
        _, outgoing, _ = massive_boundary(daughter[1:], daughter_mass_gev, incoming)
        # All validation above precedes this local mutation.
        derivative = None
        if self.cubic_history:
            candidate = self.y.copy()
            candidate[start:start+3] = daughter[1:]
            derivative = self.boundary_ray_derivative(endpoint,self.t,candidate,
                                                       mass=daughter_mass_gev,charge=0.)
        if endpoint == 0:
            self.b.jump(self.t-self.y[0], outgoing, derivative)
        else:
            self.a.jump(self.t+self.y[1], outgoing, derivative)
        self.y[start:start+3] = daughter[1:]
        self.masses[endpoint] = daughter_mass_gev
        self.charges[endpoint] = 0.
        record = {'kind': 'weak_decay', 'time_ps': self.t, 'endpoint': endpoint,
                  'position_cm': self.y[2+3*endpoint:5+3*endpoint].tolist(),
                  'parent_four_vector_gev': parent.tolist(),
                  'daughter_four_vector_gev': daughter.tolist(),
                  'parent_charge': parent_charge, 'daughter_charge': 0.,
                  'emitted_charge': parent_charge,
                  'products_gev': {k: v.tolist() for k, v in products.items() if k != 'daughter'},
                  'four_momentum_residual_gev': residual.tolist(),
                  'matching': 'LOCAL_MASSIVE_BOUNDARY_WITH_CAUSAL_NULL_RAY_JUMP'}
        self.emissions.append(record)
        return record

