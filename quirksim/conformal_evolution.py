"""Event-driven weak-decay pilot up to an explicit quantum matching surface.

The matching surface censors the classical trajectory; it does not annihilate
or reflect the pair. Nonperturbative radiation and quantum survival stay
UNKNOWN unless an independently matched quantum kernel is provided.
"""
import numpy as np
import time
from .conformal_characteristics import ConformalString
from .dynamics import sample_plan, HBARC_GEV_CM
from .physics import benchmark
from .decays import sample_decay


def evolve(representation='triplet', mass_gev=225., lambda_ev=774., seed=25,
           duration_ps=1000., mode='recoil', max_step_ps=.1,
           max_ray_angle=.003125, causal_fraction=.00625,
           matching_radius_in_confinement_units=8., observation_step_ps=5.,
           progress=None, max_steps_per_advance=200000,
           closure_stop_in_reference_units=None, max_rejected_ray_trials=None,
           minimum_step_ulps=0, wall_time_limit_seconds=None, checkpoint_path=None):
    if (mode not in ('stable', 'recoil') or duration_ps <= 0
            or not np.isfinite(duration_ps) or observation_step_ps <= 0
            or not np.isfinite(matching_radius_in_confinement_units)
            or matching_radius_in_confinement_units <= 0):
        raise ValueError('invalid evolution inputs')
    if closure_stop_in_reference_units is not None and (not np.isfinite(closure_stop_in_reference_units) or closure_stop_in_reference_units <= 0):
        raise ValueError('positive closure stop required')
    if wall_time_limit_seconds is not None and (not np.isfinite(wall_time_limit_seconds) or wall_time_limit_seconds <= 0):
        raise ValueError('positive diagnostic wall-time limit required')
    deadline = None if wall_time_limit_seconds is None else time.monotonic()+wall_time_limit_seconds
    b = benchmark(representation, mass_gev)
    solver = ConformalString.boosted_radial(b['charged_mass_gev'], lambda_ev,
                                           track_photons=True)
    initial = solver.snapshot()
    plan = sample_plan(b['lifetime_ps'], b['widths_gev'], seed)
    threshold = np.array([item['proper_time_ps'] for item in plan])
    radius = matching_radius_in_confinement_units*HBARC_GEV_CM/(lambda_ev*1e-9)
    records = [initial]
    status = 'TARGET_REACHED'
    while solver.t < duration_ps-1e-10:
        if mode == 'recoil':
            remaining = threshold-solver.y[14:16]
            for i in (0, 1):
                if solver.charges[i] != 0 and remaining[i] <= 1e-9:
                    products = sample_decay(plan[i]['channel'], solver.y[8+3*i:11+3*i],
                                            solver.masses[i], mass_gev,
                                            np.random.default_rng(plan[i]['rest_seed']))
                    event = solver.apply_endpoint_decay(i, mass_gev, products)
                    event['channel'] = plan[i]['channel']
                    event['proper_clock_residual_ps'] = float(remaining[i])
            active = [i for i in (0, 1) if solver.charges[i] != 0]
            # d(tau)/dt<=1: this cap cannot skip a proper-time threshold.
            delta = min([observation_step_ps]+[float(threshold[i]-solver.y[14+i]) for i in active])
        else:
            delta = observation_step_ps
        target = min(duration_ps, solver.t+delta)
        status = solver.advance(target, max_step_ps, causal_fraction=causal_fraction,
                                 max_ray_angle=max_ray_angle, matching_radius_cm=radius,
                                 max_steps=max_steps_per_advance,
                                 max_rejected_ray_trials=max_rejected_ray_trials,
                                 minimum_step_ulps=minimum_step_ulps, deadline_monotonic=deadline)
        records.append(solver.snapshot(samples=9))
        if (status == 'TARGET_REACHED' and closure_stop_in_reference_units is not None
                and records[-1]['chord_closure_error_cm'] >=
                closure_stop_in_reference_units*HBARC_GEV_CM/(lambda_ev*1e-9)):
            status = 'STOP_SAMPLED_CLOSURE_TOLERANCE'
            solver.last_advance_diagnostic = {'reason': status, 'time_ps': solver.t,
                'meaning': 'first failed saved observation, not exact crossing time',
                'chord_closure_error_cm': records[-1]['chord_closure_error_cm']}
        if progress is not None:
            progress({'time_ps': solver.t, 'decays': len(solver.emissions),
                      'reflections': len(solver.reflections), 'status': status,
                      'closure_cm': records[-1]['chord_closure_error_cm'],
                      'observations': len(records), 'ray_samples': len(solver.a.z)})
        if status != 'TARGET_REACHED':
            break
        if len(records) > 20000:
            status = 'STOP_OBSERVATION_BUDGET'
            break
    if checkpoint_path is not None:
        solver.save_checkpoint(checkpoint_path)
    final = solver.snapshot()
    emitted = sum((np.asarray(p, float) for e in solver.emissions
                   for p in e['products_gev'].values()), np.zeros(4))
    energy_residual = final['total_energy_gev']+emitted[0]-initial['total_energy_gev']
    momentum_residual = (np.array(final['total_momentum_gev'])+emitted[1:]
                         -final['external_impulse_gev']-initial['total_momentum_gev'])
    return {'representation': representation, 'mass_neutral_gev': mass_gev,
            'lambda_tension_ev': lambda_ev, 'seed': seed, 'mode': mode,
            'duration_requested_ps': duration_ps, 'termination': status,
            'integration': {'max_step_ps': max_step_ps, 'max_ray_angle': max_ray_angle,
                            'causal_fraction': causal_fraction},
            'last_advance_diagnostic': solver.last_advance_diagnostic,
            'checkpoint_path': str(checkpoint_path) if checkpoint_path is not None else None,
            'max_steps_per_advance': max_steps_per_advance,
            'closure_stop_in_reference_units': closure_stop_in_reference_units,
            'matching_radius_cm': radius, 'weak_plan': plan,
            'initial': initial, 'final': final, 'decays': solver.emissions,
            'wave_reflections': solver.reflections, 'records': records,
            'energy_with_weak_products_residual_gev': float(energy_residual),
            'momentum_with_weak_products_minus_field_residual_gev': momentum_residual.tolist(),
            'photon_status': 'SMOOTH_INCOHERENT_DIAGNOSTIC_NO_BACKREACTION',
            'glueball_emission_status': 'UNKNOWN_NONPERTURBATIVE_KERNEL',
            'annihilation_status': 'UNKNOWN_QUANTUM_STATE_AND_TOTAL_KERNEL',
            'conditional_on': 'NO_UNMODELED_RADIATION_OR_ANNIHILATION_BEFORE_CENSORING',
            'physical_survival_probability': None}


def quantum_handoff(classical, quantum_card):
    """Continue a censored trajectory ONLY with an explicit quantum rate card.

    The card supplies a CM level model and initial population. Subsequent
    evolution assumes a fixed pair boost and isotropic mean radiation in its
    rest frame. Recoil fluctuations and external-field evolution are absent.
    A successful calculation does not certify its physical matching.
    """
    from .radiation_annihilation import evolve_quantum_levels
    if classical['termination'] != 'QUANTUM_MATCHING_SURFACE':
        raise ValueError('a quantum matching surface is required for handoff')
    final = classical['final']
    energy = float(final['total_energy_gev'])
    momentum = np.asarray(final['total_momentum_gev'], float)
    if energy <= np.linalg.norm(momentum):
        raise ValueError('a future timelike pair four-momentum is required')
    invariant_mass = np.sqrt(energy*energy-momentum@momentum)
    levels = np.asarray(quantum_card['energies_gev'], float)
    populations = np.asarray(quantum_card['initial_probabilities'], float)
    if abs(levels@populations-invariant_mass) > 1e-8*invariant_mass:
        raise ValueError('initial quantum energy does not match the classical invariant mass')
    remaining_lab = classical['duration_requested_ps']-final['time_ps']
    if remaining_lab <= 0:
        raise ValueError('no observation time remains after matching')
    gamma = energy/invariant_mass
    result = evolve_quantum_levels(levels, quantum_card['transition_rates_per_ps'],
        quantum_card['annihilation_rates_per_ps'], quantum_card['weak_rates_per_ps'],
        populations, remaining_lab/gamma, quantum_card['rate_source'])
    mean_rad = result['expected_radiated_energy_gev']
    return {'status': 'CONDITIONAL_CLASSICAL_TO_QUANTUM_HANDOFF',
            'matching_time_ps': final['time_ps'], 'matched_invariant_mass_gev': float(invariant_mass),
            'fixed_pair_boost_gamma': float(gamma), 'quantum_evolution': result,
            'mean_radiated_lab_four_momentum_gev':
                np.r_[gamma*mean_rad, momentum/invariant_mass*mean_rad].tolist(),
            'assumptions': ['supplied quantum population and rate matrix',
                            'fixed pair boost after matching', 'isotropic CM mean radiation',
                            'weak sink ends the modeled flavor sector'],
            'physical_prediction_status': 'UNKNOWN_UNTIL_QUANTUM_CARD_IS_MATCHED'}

