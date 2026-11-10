"""Independent analytic radial/rotating references and characteristic refinement.

Run: PYTHONPATH=. python scripts/validate_conformal_boundary.py
Output is scoped engineering/analytic evidence, not a full-worldsheet PASS.
"""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import scipy
from quirksim.conformal_characteristics import ConformalString
from quirksim.conformal_radial import signed_radial_boundary
from quirksim.radial_ng_reference import ExactRadialNG
from quirksim.physics import HBAR_GEV_PS
from quirksim.dynamics import C_CM_PER_PS as C


def pulse_string(samples):
    s = np.linspace(0., 6., samples)
    angle = np.where((s > 2.) & (s < 3.), .03*np.sin(np.pi*(s-2))**4, 0.)
    a = np.c_[np.cos(angle), np.sin(angle), np.zeros_like(s)]
    b = np.tile([-1., 0, 0], (samples, 1))
    return ConformalString([225., 225.], 774., s, a, b, np.zeros((2, 3)))


def rotating_string(samples=81, beta=.2, mass=225., lambda_ev=774.):
    """Exact nonlinear rotating solution, fixed energy labels, massive ends.

    x/c=sin(omega*s)/omega*(cos(omega*t),sin(omega*t),0).
    Boundary condition: mass*gamma*beta*omega=force*cos(omega*ell).
    """
    force = (lambda_ev*1e-9)**2/HBAR_GEV_PS
    omega = force*(1-beta**2)/(mass*beta)
    ell = np.arcsin(beta)/omega
    s = np.linspace(-ell, ell, samples)
    a = np.c_[np.cos(omega*s), np.sin(omega*s), np.zeros_like(s)]
    b = np.c_[-np.cos(omega*s), np.sin(omega*s), np.zeros_like(s)]
    momentum = mass*beta/np.sqrt(1-beta**2)
    radius = C*beta/omega
    solver = ConformalString([mass, mass], lambda_ev, s, a, b,
                             [[0., -momentum, 0.], [0., momentum, 0.]],
                             [-radius, 0., 0.])
    return solver, omega, radius, momentum


def radial_check(step):
    ref = ExactRadialNG()
    segments, crossings = signed_radial_boundary(ref.mass, 774., .1, .4,
                                                1000., max_step_ps=step)
    max_position = max_momentum = max_proper = max_energy = 0.
    for sol in segments:
        for t in np.linspace(sol.t[0], sol.t[-1], 101):
            r, p, tau = sol.sol(t)
            exact = ref.state(float(t))
            max_position = max(max_position, abs(r-exact['right_position_cm']))
            max_momentum = max(max_momentum, abs(p-exact['right_momentum_gev']))
            max_proper = max(max_proper, abs(tau-exact['proper_time_each_ps']))
            energy = 2*np.hypot(ref.mass, p)+2*ref.tension*abs(r)
            max_energy = max(max_energy, abs(energy-ref.energy_total))
    expected = [ref.first_crossing_ps+i*ref.crossing_interval_ps for i in range(4)]
    crossing_error = max(abs(c['time_ps']-t) for c, t in zip(crossings, expected))
    passed = (len(crossings) == 4 and max_position < 1e-9 and max_momentum < 1e-8
              and max_proper < 1e-8 and max_energy < 1e-8 and crossing_error < 1e-8)
    return {'max_step_ps': step, 'crossings': crossings,
            'max_position_error_cm': float(max_position),
            'max_momentum_error_gev': float(max_momentum),
            'max_proper_time_error_ps': float(max_proper),
            'max_energy_error_gev': float(max_energy),
            'max_crossing_time_error_ps': float(crossing_error),
            'status': 'PASS_REDUCED_SIGNED_CONVENTION' if passed else 'FAIL'}


def pulse_check(samples, step):
    solver = pulse_string(samples)
    initial = solver.snapshot()
    status = solver.advance(9., step)
    z = solver.snapshot()
    return {'initial_samples': samples, 'max_step_ps': step,
            'integration_status': status,
            'endpoint_momenta_gev': z['momenta_gev'],
            'energy_error_gev': abs(z['total_energy_gev']-initial['total_energy_gev']),
            'momentum_error_gev': float(np.linalg.norm(np.array(z['total_momentum_gev'])
                                                       -initial['total_momentum_gev'])),
            'chord_closure_error_cm': z['chord_closure_error_cm'],
            'norm_residual': z['max_conformal_norm_residual'],
            'orthogonality_residual': z['max_conformal_orthogonality_residual']}


def rotation_check(step):
    solver, omega, radius, momentum = rotating_string()
    initial = solver.snapshot()
    errors = []
    for t in (10., 20., 40.):
        assert solver.advance(t, step) == 'TARGET_REACHED'
        z = solver.snapshot()
        angle = omega*t
        r = radius*np.array([np.cos(angle), np.sin(angle), 0.])
        p = momentum*np.array([-np.sin(angle), np.cos(angle), 0.])
        errors.append({'time_ps': t,
            'position_error_cm': float(np.max(np.abs(np.array(z['positions_cm'])-[-r, r]))),
            'momentum_error_gev': float(np.max(np.abs(np.array(z['momenta_gev'])-[-p, p]))),
            'energy_error_gev': abs(z['total_energy_gev']-initial['total_energy_gev']),
            'proper_time_error_ps': float(np.max(np.abs(np.array(z['proper_times_ps'])
                                                        -t*np.sqrt(.96)))),
            'chord_closure_error_cm': z['chord_closure_error_cm']})
    return {'max_step_ps': step, 'errors': errors}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='results/conformal_boundary_validation.json')
    args = parser.parse_args()
    report = {'scope': 'zero-field massive-boundary and characteristic prototype',
              'software': {'python': sys.version.split()[0], 'numpy': np.__version__,
                           'scipy': scipy.__version__},
              'radial': [radial_check(h) for h in (2., 1., .5)],
              'rotation': [rotation_check(h) for h in (.2, .1, .05)],
              'pulse_time_refinement': [pulse_check(241, h) for h in (.04, .02, .01)],
              'pulse_initial_grid_refinement': [pulse_check(n, .005) for n in (61, 121, 241)],
              'full_3d_continuation': 'UNKNOWN', 'physical_recoil': 'NOT_IMPLEMENTED',
              'radiation_annihilation_survival': 'UNKNOWN'}
    pulse = report['pulse_time_refinement']
    errors = [p['momentum_error_gev'] for p in pulse]
    report['pulse_time_refinement_pass'] = bool(errors[2] < errors[1] < errors[0]
                                               and errors[2] < 5e-6)
    spatial = [np.array(p['endpoint_momenta_gev'])
               for p in report['pulse_initial_grid_refinement']]
    differences = [float(np.max(np.abs(spatial[i+1]-spatial[i]))) for i in (0, 1)]
    report['pulse_grid_momentum_differences_gev'] = differences
    report['pulse_initial_grid_refinement_pass'] = differences[1] < differences[0] and differences[1] < 2e-6
    report['rotation_pass'] = all(e['position_error_cm'] < 1e-7
        and e['momentum_error_gev'] < 1e-5 and e['energy_error_gev'] < 1e-5
        for e in report['rotation'][-1]['errors'])
    report['scoped_gate'] = ('PASS' if report['rotation_pass']
        and report['pulse_time_refinement_pass']
        and report['pulse_initial_grid_refinement_pass']
        and all(r['status'] != 'FAIL' for r in report['radial']) else 'FAIL')
    path = Path(args.output); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)
    raise SystemExit(0 if report['scoped_gate'] == 'PASS' else 1)


if __name__ == '__main__':
    main()
