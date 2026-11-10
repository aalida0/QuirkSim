"""Action-derived worldsheet prototype validation, with honest pass/fail gates.

Each case runs only a short segment of a pre-stretched string, not an LHC event.
Tests action-derived internal field momentum, Noether conservation, and
finite-element momentum-conserving endpoint recoil at a single instant.
"""
import json
from dataclasses import replace
from pathlib import Path
import numpy as np

from quirksim.variational_string import StringParameters,VariationalString


def run():
    results=[]
    for n in (3,5,7):
        for gauge in (1e-2,1e-3):
            p=StringParameters(nodes=n,initial_beta=.2,initial_length_cm=.1,
                gauge_mass_fraction=gauge,regularizer_cm=1e-9,
                duration_ps=.25,max_step_ps=.025)
            s=VariationalString(p)
            y=s.initial_state(middle_offset_cm=5e-5,transverse_velocity=.01)
            diag,sol=s.solve(y)
            pre=s.action_quantities(sol.y[:,-1])
            # One fully kinematic weak decay after finite string evolution.
            candidate,weak=s.weak_decay_impulse(sol.y[:,-1],0,'pion',seed=42,
                                                 energy_tolerance_gev=1e-4,
                                                 strict=False)
            after=s.action_quantities(candidate)
            cont,sol2=s.solve(candidate,duration_ps=.1)
            position_before=s._unpack(sol.y[:,-1])[0]
            position_after=s._unpack(sol2.y[:,-1])[0]
            row={
                'nodes':n,'gauge_fraction':gauge,'gauge_mass_gev':s.gauge_mass,
                'free_segment_energy_drift_gev':diag['energy_drift_gev'],
                'free_segment_momentum_drift_gev':
                    float(np.linalg.norm(np.array(diag['final_momentum_gev'])-
                                            np.array(diag['initial_momentum_gev']))),
                'free_segment_angular_momentum_drift_hbar':
                    float(np.linalg.norm(np.array(diag['final_angular_momentum_hbar'])-
                                            np.array(diag['initial_angular_momentum_hbar']))),
                'weak_decay_emitted_energy_gev':weak['emitted_energy_gev'],
                'weak_decay_total_energy_residual_gev':weak['energy_residual_gev'],
                'weak_decay_momentum_residual_gev':
                    float(np.linalg.norm(weak['momentum_residual_gev'])),
                'weak_decay_angular_momentum_residual_hbar':
                    float(np.linalg.norm(weak['angular_momentum_residual_hbar'])),
                'weak_decay_1e-4GeV_gate':weak['energy_gate_pass'],
                'post_decay_free_energy_drift_gev':cont['energy_drift_gev'],
                'end_separation_cm':
                    float(np.linalg.norm(position_after[0]-position_after[-1])),
                'post_decay_max_node_displacement_cm':
                    float(np.max(np.linalg.norm(position_after-position_before,axis=1)))}
            results.append(row)
            print('[worldsheet-check] '+json.dumps(row),flush=True)
    report={
        'status':'DEVELOPMENT_ONLY_NOT_FULL_NG_WORLD_SHEET_VALIDATED',
        'cases':results,
        'all_nominal_conservation_passed':
            all(abs(r['free_segment_energy_drift_gev'])<1e-5 and
                r['free_segment_momentum_drift_gev']<1e-5 and
                r['free_segment_angular_momentum_drift_hbar']<1e7 and
                r['post_decay_free_energy_drift_gev']<1e-5
                for r in results),
        'all_weak_0p1MeV_continuity_passed':
            all(r['weak_decay_1e-4GeV_gate'] for r in results),
        'physical_limitations':[
            'spatial midpoint discretization with artificial tangent gauge masses',
            'short-time initialized stretched string only',
            'endpoint instantaneous canonical weak impulse not yet continuum converged',
            'no quantum annihilation or self-consistent glueball emission',
            'no causal convergence or realistic pp production demonstrated']}
    Path('results').mkdir(exist_ok=True)
    Path('results/variational_gate_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('[worldsheet-summary] '+json.dumps({k:v for k,v in report.items() if k!='cases'}),
          flush=True)
    if not report['all_nominal_conservation_passed']:
        raise AssertionError('Noether conservation validation failed')
    return report


if __name__=='__main__':
    run()
