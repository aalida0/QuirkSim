"""Causal magnetic/weak pilot with explicit short-distance censoring."""
import argparse
import json
from pathlib import Path
from quirksim.conformal_evolution import evolve
from quirksim.dynamics import HBARC_GEV_CM


def validate(result):
    ref = HBARC_GEV_CM/(result['lambda_tension_ev']*1e-9)
    closure = max(r['chord_closure_error_cm'] for r in result['records'])
    energy = abs(result['energy_with_weak_products_residual_gev'])
    momentum = max(abs(x) for x in result['momentum_with_weak_products_minus_field_residual_gev'])
    passed = (result['termination'] in ('TARGET_REACHED', 'QUANTUM_MATCHING_SURFACE')
              and closure < .1*ref and energy < 1e-7 and momentum < 5e-6)
    result['numerical_validation'] = {'status': 'PASS' if passed else 'FAIL',
        'max_chord_closure_error_cm': closure,
        'max_closure_over_confinement_reference': closure/ref,
        'energy_tolerance_gev': 1e-7, 'momentum_tolerance_gev': 5e-6,
        'closure_tolerance_in_confinement_units': .1,
        'meaning': 'SELECTED_CLASSICAL_EVOLUTION_ONLY_NOT_QUANTUM_SURVIVAL'}
    return passed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mode', choices=['stable', 'recoil'], default='recoil')
    p.add_argument('--model', choices=['triplet', 'doublet'], default='triplet')
    p.add_argument('--mass', type=float)
    p.add_argument('--lambda-ev', type=float)
    p.add_argument('--seed', type=int, default=25)
    p.add_argument('--duration-ps', type=float, default=1000.)
    p.add_argument('--max-step-ps', type=float, default=.05)
    p.add_argument('--ray-angle', type=float, default=.00078125)
    p.add_argument('--causal-fraction', type=float, default=.0015625)
    p.add_argument('--max-refinements', type=int, choices=range(4), default=1)
    p.add_argument('--output', default='results/conformal_evolution.json')
    a = p.parse_args()
    last_print = [-50.]
    def progress(row):
        if row['time_ps'] >= last_print[0]+50. or row['status'] != 'TARGET_REACHED':
            print(json.dumps(row), flush=True)
            last_print[0] = row['time_ps']
    path = Path(a.output); path.parent.mkdir(parents=True, exist_ok=True)
    attempts = []
    for refinement in range(a.max_refinements+1):
        last_print[0] = -50.
        factor = 2**refinement
        r = evolve(representation=a.model,
                   mass_gev=a.mass if a.mass is not None else (225. if a.model == 'triplet' else 150.),
                   lambda_ev=a.lambda_ev if a.lambda_ev is not None else (774. if a.model == 'triplet' else 2154.),
                   mode=a.mode, seed=a.seed, duration_ps=a.duration_ps,
                   max_step_ps=a.max_step_ps/factor, max_ray_angle=a.ray_angle/factor,
                   causal_fraction=a.causal_fraction/factor, progress=progress)
        passed = validate(r)
        attempts.append({'refinement': refinement, 'integration': r['integration'],
                         'validation': r['numerical_validation'], 'termination': r['termination']})
        r['refinement_attempts'] = attempts
        path.write_text(json.dumps(r, indent=2)+'\n')
        if passed or r['termination'] not in ('TARGET_REACHED', 'QUANTUM_MATCHING_SURFACE'):
            break
    print(json.dumps({k: r[k] for k in ['mode', 'seed', 'termination',
                    'energy_with_weak_products_residual_gev', 'numerical_validation']}), flush=True)
    raise SystemExit(0 if passed else 1)


if __name__ == '__main__':
    main()
