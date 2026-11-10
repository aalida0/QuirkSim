"""Causal magnetic/weak pilot with explicit short-distance censoring."""
import argparse
import json
from pathlib import Path
import numpy as np
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
    p.add_argument('--seed', type=int, default=25)
    p.add_argument('--duration-ps', type=float, default=1000.)
    p.add_argument('--max-step-ps', type=float, default=.05)
    p.add_argument('--ray-angle', type=float, default=.00078125)
    p.add_argument('--causal-fraction', type=float, default=.0015625)
    p.add_argument('--output', default='results/conformal_evolution.json')
    a = p.parse_args()
    r = evolve(mode=a.mode, seed=a.seed, duration_ps=a.duration_ps,
               max_step_ps=a.max_step_ps, max_ray_angle=a.ray_angle,
               causal_fraction=a.causal_fraction)
    passed = validate(r)
    path = Path(a.output); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(r, indent=2)+'\n')
    print(json.dumps({k: r[k] for k in ['mode', 'seed', 'termination',
                    'energy_with_weak_products_residual_gev', 'numerical_validation']}), flush=True)
    raise SystemExit(0 if passed else 1)


if __name__ == '__main__':
    main()
