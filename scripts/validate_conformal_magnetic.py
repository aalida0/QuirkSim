"""Reproduce the independent 2 T FEM comparison and matching-surface refinement."""
import argparse
import json
from pathlib import Path
import numpy as np
from quirksim.conformal_characteristics import ConformalString
from quirksim.ng_quadrature import GaussQuadratureNGString
from quirksim.variational_string import StringParameters
from quirksim.physics import benchmark
from quirksim.dynamics import HBARC_GEV_CM

REFERENCE_CM = HBARC_GEV_CM/(774e-9)


def compare_fem():
    rows = []
    for nodes in (7, 11, 15):
        fem = GaussQuadratureNGString(StringParameters(nodes=nodes, b_tesla=2.,
            max_step_ps=.5, regularizer_cm=2e-10), 2)
        initial = fem.boosted_straight_initial_state(.25)
        conformal = ConformalString.boosted_radial(fem.masses[0], 774.)
        for t in (20., 100., 200.):
            _, sol = fem.solve(initial, duration_ps=t)
            positions, _ = fem._unpack(sol.y[:, -1])
            status = conformal.advance(t, .05)
            if status != 'TARGET_REACHED':
                raise RuntimeError(status)
            snap = conformal.snapshot()
            difference = float(np.max(np.abs(np.array(snap['positions_cm'])-positions[[0, -1]])))
            rows.append({'nodes': nodes, 'time_ps': t,
                'max_position_difference_cm': difference,
                'fem_positions_cm': positions[[0, -1]].tolist(),
                'conformal_positions_cm': snap['positions_cm']})
            print(json.dumps(rows[-1]), flush=True)
    passed = all(next(r['max_position_difference_cm'] for r in rows if r['nodes']==15 and r['time_ps']==t)
                 < next(r['max_position_difference_cm'] for r in rows if r['nodes']==11 and r['time_ps']==t)
                 < next(r['max_position_difference_cm'] for r in rows if r['nodes']==7 and r['time_ps']==t)
                 for t in (20., 100., 200.))
    passed &= max(r['max_position_difference_cm'] for r in rows if r['nodes']==15) < .1*REFERENCE_CM
    return {'records': rows, 'status': 'PASS' if passed else 'FAIL',
            'scope': 'selected 2 T longitudinally boosted stretched triplet at common lab times'}


def matching_refinement():
    rows = []
    for fraction, angle in ((.00625, .003125), (.003125, .0015625), (.0015625, .00078125)):
        s = ConformalString.boosted_radial(benchmark()['charged_mass_gev'], 774.)
        initial = s.snapshot()
        status = s.advance(1000., .05, minimum_interval_ps=1e-12,
            causal_fraction=fraction, max_ray_angle=angle, matching_radius_cm=8*REFERENCE_CM)
        snap = s.snapshot()
        row = {'causal_fraction': fraction, 'max_ray_angle': angle, 'termination': status,
               'initial': initial, 'final': snap,
               'poincare_129': s.poincare(129), 'poincare_257': s.poincare(257)}
        rows.append(row)
        print(json.dumps({'fraction': fraction, 'termination': status,
                          'time_ps': s.t, 'closure_cm': snap['chord_closure_error_cm']}), flush=True)
    errors = [r['final']['chord_closure_error_cm'] for r in rows]
    passed = (all(r['termination']=='QUANTUM_MATCHING_SURFACE' for r in rows)
              and errors[2] < errors[1] < errors[0] and errors[2] < .1*REFERENCE_CM
              and abs(rows[-1]['final']['time_ps']-rows[-2]['final']['time_ps']) < 5e-7)
    return {'records': rows, 'status': 'PASS' if passed else 'FAIL',
            'scope': 'matching radius 8 hbar*c/Lambda; not quantum survival or contact continuation'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--suite', choices=['fem', 'matching', 'both'], default='both')
    p.add_argument('--output', default='results/conformal_magnetic_validation.json')
    a = p.parse_args()
    report = {'scope': 'classical magnetic cross-validation only'}
    if a.suite in ('fem', 'both'):
        report['fem_comparison'] = compare_fem()
    if a.suite in ('matching', 'both'):
        report['matching_refinement'] = matching_refinement()
    report['status'] = 'PASS' if all(v['status']=='PASS' for v in report.values() if isinstance(v, dict)) else 'FAIL'
    path = Path(a.output); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2)+'\n')
    raise SystemExit(0 if report['status']=='PASS' else 1)


if __name__ == '__main__':
    main()
