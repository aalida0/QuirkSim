"""Scoped radiation/annihilation audit and explicit assumption study.

The output is never a physical survival prediction or a calibrated glueball
rate. The alpha_IC value, overlap and emission efficiencies are supplied
nuisance assumptions; they are not extracted from the classical geometry.
"""
import argparse
import json
import math
from pathlib import Path
from quirksim.physics import benchmark
from quirksim.dynamics import HBARC_GEV_CM
from quirksim.radiation_annihilation import (
    OrbitAveragedModel, dirac_singlet_partial_sigma_v, evolve_orbit_average,
    linear_cooling_time_ps, radiation_recoil_scale, coulomb_ground_reference,
    evolve_quantum_levels)


def audit():
    cards, scenarios = [], []
    for representation, mass, lam in [('triplet', 225., 774.), ('doublet', 150., 2154.)]:
        b = benchmark(representation, mass)
        m = b['charged_mass_gev']
        excess = 2*m*(1/math.sqrt(1-.4**2)-1)+(lam*1e-9)**2/HBARC_GEV_CM*.1
        for n in (2, 3):
            kernel = dirac_singlet_partial_sigma_v(m, n, .04)
            card = {'representation': representation, 'benchmark': b,
                    'lambda_tension_ev': lam, 'n_ic_assumed': n,
                    'alpha_ic_assumed': .04, 'initial_pair_excess_energy_gev': excess,
                    'initial_rest_length_cm': .1, 'initial_beta': .4,
                    'partial_kernel': kernel}
            # An explicitly supplied Coulomb potential is an analytic reference,
            # not a claim that the hard-scale coupling matches the Bohr scale.
            strength = (n*n-1)/(2*n)*.04
            ground = coulomb_ground_reference(m, m, strength)
            width = ground['psi_squared_gev3']*kernel['sum_partial_gev_minus2']
            from quirksim.physics import HBAR_GEV_PS
            card['ground_state_reference'] = {**ground,
                'potential_strength_assumed': strength, 'partial_ensemble_average_width_gev': width,
                'partial_ensemble_average_lifetime_ps': HBAR_GEV_PS/width,
                'spin_resolved_physical_lifetime': None}
            cards.append(card)
            for epsilon in (0., .1, 1.):
                for overlap in (0., .1, 1.):
                    model = OrbitAveragedModel(m, lam, epsilon, lam*1e-9,
                        kernel['sum_partial_gev_minus2'], overlap, b['lifetime_ps'],
                        'Li2024 Table II + explicit photon-to-fermion thresholds; alpha_IC=0.04')
                    result = evolve_orbit_average(model, excess, 1000., 8*lam*1e-9)
                    cooling = linear_cooling_time_ps(model, excess, 8*lam*1e-9)
                    result.update({'representation': representation, 'n_ic_assumed': n,
                        'linear_cooling_time_ps_conditional': cooling,
                        'cooling_time_over_weak_lifetime': cooling/b['lifetime_ps'],
                        'radiation_recoil_materiality': radiation_recoil_scale(
                            result['conditional_alive_trajectory_energy_loss_gev'], m, 1000., lam)})
                    scenarios.append(result)
    # The supplied test chain verifies a complete radiation->annihilation
    # calculation. It is deliberately not assigned to the physical benchmarks.
    chain = evolve_quantum_levels([451., 450.], [[0., 2.], [0., 0.]],
                [0., 1.], [0., 0.], [1., 0.], 4., 'analytic validation chain ONLY')
    return {'implementation_status': 'CONDITIONAL_RADIATION_AND_ANNIHILATION_LAYERS_IMPLEMENTED',
            'physical_prediction_status': 'UNKNOWN', 'benchmark_cards': cards,
            'conditional_orbit_averaged_scenarios': scenarios,
            'analytic_quantum_cascade_validation': chain,
            'conventions': ['Lambda_tension=sqrt(string tension), not an inferred Lambda_MS',
                'one full signed period=two radial crossings',
                'glueball loss per full signed period=epsilon times supplied emission energy',
                'orbit-averaged rates concern the state before its first weak decay',
                'the example emission energy Lambda is a phenomenological energy-loss prescription, not an on-shell glueball spectrum'],
            'mandatory_missing_inputs': ['nonperturbative glueball transition and recoil kernel',
                'quantum spin/flavor/partial-wave state at the matching surface',
                'complete EW short-distance annihilation matrix',
                'matched quantum-level transition and capture rates'],
            'release': 'NO_PHYSICAL_SURVIVAL_OR_RECOIL_SUPPRESSION_CLAIM'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='results/radiation_annihilation_audit.json')
    args = parser.parse_args()
    report = audit()
    path = Path(args.output); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'output': str(path), 'status': report['implementation_status'],
                      'physical_prediction_status': report['physical_prediction_status'],
                      'scenarios': len(report['conditional_orbit_averaged_scenarios'])}), flush=True)


if __name__ == '__main__':
    main()
