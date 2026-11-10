"""Event-driven weak-decay pilot up to an explicit quantum matching surface.

The matching surface censors the classical trajectory; it does not annihilate
or reflect the pair. Nonperturbative radiation and quantum survival stay
UNKNOWN unless an independently matched quantum kernel is provided.
"""
import numpy as np
from .conformal_characteristics import ConformalString
from .dynamics import sample_plan, HBARC_GEV_CM
from .physics import benchmark
from .decays import sample_decay


def evolve(representation='triplet', mass_gev=225., lambda_ev=774., seed=25,
           duration_ps=1000., mode='recoil', max_step_ps=.1,
           max_ray_angle=.003125, causal_fraction=.00625,
           matching_radius_in_confinement_units=8., observation_step_ps=5.):
    if (mode not in ('stable', 'recoil') or duration_ps <= 0
            or not np.isfinite(duration_ps) or observation_step_ps <= 0
            or not np.isfinite(matching_radius_in_confinement_units)
            or matching_radius_in_confinement_units <= 0):
        raise ValueError('invalid evolution inputs')
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
                                 max_ray_angle=max_ray_angle, matching_radius_cm=radius)
        records.append(solver.snapshot(samples=9))
        if status != 'TARGET_REACHED':
            break
        if len(records) > 20000:
            status = 'STOP_OBSERVATION_BUDGET'
            break
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
