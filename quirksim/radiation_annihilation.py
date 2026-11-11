"""Explicitly conditional radiation and short-distance annihilation machinery.

No classical radius determines a quantum s-wave occupation. No calibrated
nonperturbative glueball kernel is supplied. This module provides traceable
short-distance partial coefficients, smooth photon power, and a conservative
orbit-averaged competing-risk calculation with declared inputs. It never
silently changes a moving worldsheet to remove energy.
"""
from dataclasses import dataclass, asdict
import math
import numpy as np
from scipy.integrate import quad, solve_ivp
from .physics import HBAR_GEV_PS
from .dynamics import C_CM_PER_PS, HBARC_GEV_CM

ALPHA_ZERO = 1/137.035999084


def lienard_power_gev_per_ps(beta, beta_dot_per_ps, charge=1., alpha=ALPHA_ZERO):
    """Single smooth charge, rationalized EM: P=(2 alpha/3) hbar gamma^6(...).

    This excludes impulsive weak-vertex bremsstrahlung and interference between
    the endpoints. It is a radiation diagnostic, not a radiation-reaction force.
    """
    v, a = np.asarray(beta, float), np.asarray(beta_dot_per_ps, float)
    if (v.shape != (3,) or a.shape != (3,) or not np.all(np.isfinite(v))
            or not np.all(np.isfinite(a)) or not math.isfinite(charge)
            or not math.isfinite(alpha) or alpha <= 0 or v@v >= 1):
        raise ValueError('finite subluminal velocity and physical EM inputs required')
    g2 = 1/(1-float(v@v))
    # Positive decomposition avoids cancellation for ultrarelativistic motion.
    invariant = (1-float(v@v))*float(a@a)+float(v@a)**2
    return (2/3)*alpha*charge**2*HBAR_GEV_PS*g2**3*invariant


def dirac_singlet_partial_sigma_v(mass_gev, n_ic, alpha_ic, charge=1., alpha_em=ALPHA_ZERO):
    """Color-neutral Dirac pair, IC singlet, unpolarized hard threshold kernels.

    Li et al. PRD109(2024)095005 Table II: gamma-gamma, e-e and IC gluons.
    Other charged fermions via photon exchange include threshold factors.
    Weak bosons, Z exchange/interference, Sommerfeld and off-diagonal flavor
    channels are NOT included. These coefficients are not an EW-multiplet total.
    """
    if (not all(map(math.isfinite, [mass_gev, alpha_ic, charge, alpha_em]))
            or mass_gev <= 0 or not isinstance(n_ic, int) or n_ic < 2
            or not 0 < alpha_ic < 1 or not 0 < alpha_em < 1):
        raise ValueError('positive mass and perturbative couplings, integer NIC>=2 required')
    base = math.pi/mass_gev**2
    result = {'ic_gluons': base*alpha_ic**2*(n_ic*n_ic-1)/(4*n_ic),
              'gamma_gamma': base*n_ic*alpha_em**2*charge**4}
    # Fixed illustrative SM pole thresholds; these are not precision EW inputs.
    fermions = [('electron', .00051099895, 1., -1.), ('muon', .1056583755, 1., -1.),
                ('tau', 1.77686, 1., -1.), ('up', .0022, 3., 2/3),
                ('down', .0047, 3., -1/3), ('strange', .096, 3., -1/3),
                ('charm', 1.67, 3., 2/3), ('bottom', 4.78, 3., -1/3),
                ('top', 172.5, 3., 2/3)]
    for name, mf, colors, qf in fermions:
        r = mf/mass_gev
        result['photon_to_'+name] = (base*n_ic*alpha_em**2*charge**2*colors*qf**2
            *math.sqrt(1-r*r)*(1+r*r/2) if r < 1 else 0.)
    return {'sigma_v_gev_minus2': result, 'sum_partial_gev_minus2': sum(result.values()),
            'coverage': 'PARTIAL_QED_AND_IC_ONLY',
            'state': 'SAME_FLAVOR_DIRAC_UNPOLARIZED_IC_SINGLET',
            'source': '10.1103/PhysRevD.109.095005 Table II',
            'missing': ['EW Z/W channels and interference', 'spin/flavor state matching',
                        'Sommerfeld/coupled-channel matching', 'nonperturbative overlap']}


def wkb_crossing_exponent(mass_gev, beta, sigma_v_gev_minus2, s_wave_overlap):
    """Kang-Luty Eq.(5.29) as a small-probability optical-depth increment.

    The overlap is a supplied occupation/penetration factor, never inferred
    from a classical endpoint distance or total angular momentum.
    """
    if (not all(map(math.isfinite, [mass_gev, beta, sigma_v_gev_minus2, s_wave_overlap]))
            or mass_gev <= 0 or not 0 <= beta < 1 or sigma_v_gev_minus2 < 0
            or not 0 <= s_wave_overlap <= 1):
        raise ValueError('invalid short-distance crossing inputs')
    return s_wave_overlap*mass_gev**2*beta*sigma_v_gev_minus2/(2*math.pi)


@dataclass(frozen=True)
class OrbitAveragedModel:
    mass_gev: float
    lambda_tension_ev: float
    epsilon_per_full_signed_period: float
    energy_per_emission_gev: float
    sigma_v_gev_minus2: float
    s_wave_overlap: float
    endpoint_weak_lifetime_ps: float
    kernel_source: str
    kernel_coverage: str = 'partial'
    charge_magnitude: float = 1.
    include_photons: bool = True

    def __post_init__(self):
        numbers = [self.mass_gev, self.lambda_tension_ev,
                   self.epsilon_per_full_signed_period, self.energy_per_emission_gev,
                   self.sigma_v_gev_minus2, self.s_wave_overlap,
                   self.endpoint_weak_lifetime_ps, self.charge_magnitude]
        if (not all(map(math.isfinite, numbers)) or self.mass_gev <= 0
                or self.lambda_tension_ev <= 0 or self.energy_per_emission_gev <= 0
                or not 0 <= self.epsilon_per_full_signed_period <= 1
                or self.sigma_v_gev_minus2 < 0 or not 0 <= self.s_wave_overlap <= 1
                or self.endpoint_weak_lifetime_ps <= 0 or not self.kernel_source.strip()
                or self.kernel_coverage not in ('partial', 'provided_total')):
            raise ValueError('all physics assumptions must be finite, explicit and physical')

    def rates(self, pair_excess_energy_gev):
        m, exc = self.mass_gev, pair_excess_energy_gev
        if not math.isfinite(exc) or exc <= 0:
            raise ValueError('positive pair excess energy required in the linear regime')
        energy_each = m+exc/2
        p = math.sqrt(exc*(m+exc/4))
        beta = p/energy_each
        force = (self.lambda_tension_ev*1e-9)**2/HBAR_GEV_PS
        crossing_period = 2*p/force
        full_signed_period = 2*crossing_period
        glue_power = self.epsilon_per_full_signed_period*self.energy_per_emission_gev/full_signed_period
        photon_power = ((4/3)*ALPHA_ZERO*self.charge_magnitude**2
                        *HBAR_GEV_PS*(force/m)**2 if self.include_photons else 0.)
        optical_depth = wkb_crossing_exponent(m, beta, self.sigma_v_gev_minus2, self.s_wave_overlap)
        average_proper_rate = m/p*math.asinh(p/m)
        return {'crossing_period_ps': crossing_period, 'full_signed_period_ps': full_signed_period,
                'glueball_power_gev_per_ps': glue_power,
                'photon_power_gev_per_ps': photon_power,
                'annihilation_rate_per_ps': optical_depth/crossing_period,
                'first_weak_rate_per_ps': 2*average_proper_rate/self.endpoint_weak_lifetime_ps,
                'crossing_optical_depth': optical_depth}


def evolve_orbit_average(model, initial_excess_energy_gev, duration_ps,
                         linear_stop_energy_gev, max_step_ps=10.):
    """Radiative energy feedback and competing first-weak/annihilation risks.

    Radial zero-field *orbit average*, not a continuation of a curved string.
    Stops at the supplied linear-to-quantum matching energy. No unsupported
    ground-state capture, invisible BR, or return to a worldsheet is installed.
    """
    if (not all(map(math.isfinite, [initial_excess_energy_gev, duration_ps,
                                   linear_stop_energy_gev, max_step_ps]))
            or not 0 < linear_stop_energy_gev < initial_excess_energy_gev
            or duration_ps <= 0 or max_step_ps <= 0):
        raise ValueError('positive time/step and ordered energy matching scales required')

    def rhs(t, y):
        rates = model.rates(max(y[0], linear_stop_energy_gev))
        power = rates['glueball_power_gev_per_ps']+rates['photon_power_gev_per_ps']
        ann, weak = rates['annihilation_rate_per_ps'], rates['first_weak_rate_per_ps']
        return [-power, -(ann+weak)*y[1], ann*y[1], weak*y[1], power*y[1]]

    def stop(t, y):
        return y[0]-linear_stop_energy_gev

    stop.terminal = True
    stop.direction = -1
    sol = solve_ivp(rhs, (0., duration_ps), [initial_excess_energy_gev, 1., 0., 0., 0.],
                    rtol=2e-10, atol=1e-13, method='DOP853', max_step=max_step_ps, events=stop)
    if not sol.success:
        raise RuntimeError(sol.message)
    energy, survival, ann, weak, radiation = sol.y[:, -1]
    probabilities = [float(survival), float(ann), float(weak)]
    if min(probabilities) < -1e-10 or abs(sum(probabilities)-1) > 1e-9:
        raise RuntimeError('probability conservation failed')
    return {'model': asdict(model), 'status': 'CONDITIONAL_ORBIT_AVERAGED_MODEL_ONLY',
            'termination': 'LINEAR_TO_QUANTUM_MATCHING' if len(sol.t_events[0]) else 'TIME_LIMIT',
            'time_ps': float(sol.t[-1]), 'pair_excess_energy_gev': float(energy),
            'survival_before_first_weak_or_annihilation': float(survival),
            'annihilation_before_first_weak_probability': float(ann),
            'first_weak_before_annihilation_probability': float(weak),
            'probability_residual': float(sum(probabilities)-1),
            'conditional_alive_trajectory_energy_loss_gev': float(initial_excess_energy_gev-energy),
            'survival_weighted_expected_radiated_energy_gev': float(radiation),
            'linear_stop_energy_gev': linear_stop_energy_gev,
            'physical_prediction_status': 'UNKNOWN',
            'missing_for_physical_prediction': ['calibrated glueball emission and recoil kernel',
                'quantum spin/flavor/partial-wave density matrix at matching',
                'full short-distance EW kernel', 'nonradial and magnetic matching',
                'quantum-regime capture/cooling evolution']}


def linear_cooling_time_ps(model, initial_energy_gev, stop_energy_gev):
    if not 0 < stop_energy_gev < initial_energy_gev:
        raise ValueError('invalid cooling energy interval')
    initial_rates = model.rates(initial_energy_gev)
    if initial_rates['glueball_power_gev_per_ps']+initial_rates['photon_power_gev_per_ps'] == 0:
        return math.inf
    return quad(lambda e: 1/(model.rates(e)['glueball_power_gev_per_ps']
                              +model.rates(e)['photon_power_gev_per_ps']),
                stop_energy_gev, initial_energy_gev, epsabs=1e-6, epsrel=1e-9)[0]


def radiation_recoil_scale(total_emitted_energy_gev, endpoint_mass_gev, horizon_ps, lambda_ev):
    """A recoil materiality *scale*, not a bound on nonlinear string response.

    Assigning all emitted momentum to one heavy endpoint gives a useful free
    flight displacement scale. It cannot bound phase errors or a near-contact
    scattering trajectory without an actual local radiation kernel.
    """
    if (not all(map(math.isfinite, [total_emitted_energy_gev, endpoint_mass_gev, horizon_ps, lambda_ev]))
            or total_emitted_energy_gev < 0 or min(endpoint_mass_gev, horizon_ps, lambda_ev) <= 0):
        raise ValueError('invalid radiation scale inputs')
    displacement = C_CM_PER_PS*horizon_ps*total_emitted_energy_gev/endpoint_mass_gev
    reference = HBARC_GEV_CM/(lambda_ev*1e-9)
    return {'free_flight_recoil_displacement_scale_cm': displacement,
            'scale_over_confinement_reference': displacement/reference,
            'worldsheet_geometry_error_bound': None,
            'interpretation': 'SMALL_ENERGY_LOSS_DOES_NOT_CERTIFY_SHORT_DISTANCE_GEOMETRY'}


def coulomb_ground_reference(mass1_gev, mass2_gev, potential_strength):
    """Analytic -kappa/r Schrödinger reference; not a confinement-scale match."""
    if (not all(map(math.isfinite, [mass1_gev, mass2_gev, potential_strength]))
            or min(mass1_gev, mass2_gev, potential_strength) <= 0 or potential_strength >= 1):
        raise ValueError('positive masses and perturbative potential strength required')
    reduced_mass = mass1_gev*mass2_gev/(mass1_gev+mass2_gev)
    return {'binding_energy_gev': reduced_mass*potential_strength**2/2,
            'psi_squared_gev3': (reduced_mass*potential_strength)**3/math.pi,
            'status': 'ANALYTIC_COULOMB_REFERENCE_ONLY'}


def evolve_quantum_levels(energies_gev, transition_rates_per_ps,
                           annihilation_rates_per_ps, weak_rates_per_ps,
                           initial_probabilities, duration_ps, rate_source):
    """Trace-preserving finite-level radiation/annihilation master equation.

    Transition matrix [i,j] is the supplied rate from i to j. Only downward
    radiation is allowed. All quantum matching, spin/flavor populations and
    transition/annihilation rates must be provided with provenance. Energies
    include the pair rest energy. Radiation and absorbing channels have
    separate probability and mean-energy ledgers.
    """
    e = np.asarray(energies_gev, float)
    rates = np.asarray(transition_rates_per_ps, float)
    ann, weak = np.asarray(annihilation_rates_per_ps, float), np.asarray(weak_rates_per_ps, float)
    initial = np.asarray(initial_probabilities, float)
    n = len(e)
    if (e.ndim != 1 or n < 1 or rates.shape != (n, n) or ann.shape != (n,)
            or weak.shape != (n,) or initial.shape != (n,) or not rate_source.strip()
            or not all(np.all(np.isfinite(v)) for v in [e, rates, ann, weak, initial])
            or np.any(e <= 0) or np.any(rates < 0) or np.any(ann < 0)
            or np.any(weak < 0) or np.any(initial < 0)
            or abs(initial.sum()-1) > 1e-12 or np.any(np.diag(rates) != 0)
            or np.any(rates[(e[:, None]-e[None, :]) <= 0] != 0)
            or not math.isfinite(duration_ps) or duration_ps <= 0):
        raise ValueError('normalized state, positive physical downward rates and provenance required')
    generator = rates.T-np.diag(rates.sum(axis=1)+ann+weak)
    photon_energy_rate = (rates*(e[:, None]-e[None, :])).sum(axis=1)

    # Linear constant-rate system: matrix exponential avoids stiffness and
    # preserves the exact analytic model to numerical roundoff.
    from scipy.linalg import expm
    matrix = np.zeros((n+5, n+5))
    matrix[:n, :n] = generator
    matrix[n:, :n] = np.array([ann, weak, photon_energy_rate, ann*e, weak*e])
    final = expm(matrix*duration_ps)@np.r_[initial, np.zeros(5)]
    live, pann, pweak, radiation, eann, eweak = final[:n], *final[n:]
    probability_residual = live.sum()+pann+pweak-1
    energy_residual = e@live+radiation+eann+eweak-e@initial
    if min(np.r_[live, pann, pweak, radiation, eann, eweak]) < -1e-10:
        raise RuntimeError('nonpositive quantum evolution')
    if abs(probability_residual) > 1e-9 or abs(energy_residual) > 1e-9*max(e):
        raise RuntimeError('quantum probability/energy accounting failed')
    return {'status': 'CONDITIONAL_FINITE_LEVEL_QUANTUM_MODEL', 'rate_source': rate_source,
            'time_ps': duration_ps, 'level_probabilities': live.tolist(),
            'annihilation_probability': float(pann), 'first_weak_probability': float(pweak),
            'expected_radiated_energy_gev': float(radiation),
            'expected_energy_to_annihilation_gev': float(eann),
            'expected_energy_to_weak_sink_gev': float(eweak),
            'probability_residual': float(probability_residual),
            'energy_residual_gev': float(energy_residual),
            'physical_prediction_status': 'UNKNOWN_UNTIL_RATES_AND_INITIAL_DENSITY_ARE_MATCHED'}
