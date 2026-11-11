"""First-return / 260 ps long-time mesh-health stress test, no weak decays.

THIS IS NOT a validated physical quirk simulation. Failing or collapsed
meshes are flagged, not converted into geometric annihilation results.
"""
import json
from pathlib import Path
import numpy as np
from quirksim.variational_string import StringParameters, VariationalString
from quirksim.worldsheet_experiment import pair_boost_initial_state


def single(field_tesla):
    p=StringParameters(nodes=5,mass_gev=225,lambda_ev=774,
                       initial_length_cm=.1,initial_beta=.4,
                       b_tesla=field_tesla,max_step_ps=.5,
                       gauge_mass_fraction=1e-3,
                       duration_ps=260)
    solver=VariationalString(p)
    y=pair_boost_initial_state(solver,.25)
    try:
        diag,sol=solver.solve(y,duration_ps=260)
        nodes=p.nodes
        positions=sol.y[:3*nodes,:].T.reshape(-1,nodes,3)
        spans=np.linalg.norm(np.diff(positions,axis=1),axis=2)
        min_ell=float(spans.min())
        ratio=float(max(np.max(spans,axis=1)/
                        np.maximum(np.min(spans,axis=1),p.regularizer_cm)))
        end_distance=np.linalg.norm(positions[:,0,:]-positions[:,-1,:],axis=1)
        return {
         'status':'NUMERICAL_STRESS_RUN_COMPLETE_NOT_PHYSICS_VALIDATION',
         'field_tesla':field_tesla,
         'minimum_element_length_cm':min_ell,
         'worst_element_length_ratio':ratio,
         'mesh_collapse_under_100reg':min_ell<100*p.regularizer_cm,
         'mesh_imbalance_above_50':ratio>50,
         'min_endpoint_separation_cm':float(end_distance.min()),
         'max_endpoint_separation_cm':float(end_distance.max()),
         'initial_energy_gev':diag['initial_energy_gev'],
         'energy_drift_gev':diag['energy_drift_gev'],
         'initial_momentum_gev':diag['initial_momentum_gev'],
         'final_momentum_gev':diag['final_momentum_gev'],
         'initial_L_over_hbar':diag['initial_angular_momentum_hbar'],
         'final_L_over_hbar':diag['final_angular_momentum_hbar'],
         'nfev':diag['n_rhs_evaluations'],
         'n_accepted_time_samples':len(sol.t),
         'reliability':'INSUFFICIENT' if min_ell<100*p.regularizer_cm or ratio>50
                       else 'STILL_REQUIRES_MESH_REFINEMENT'
        }
    except Exception as exc:
        return {'status':'BLOCKED_BY_NUMERICAL_STRING_EVOLUTION',
                'field_tesla':field_tesla,
                'reason':str(exc)}


def run():
    reports=[single(0.),single(2.)]
    result={'reports':reports,
            'scientific_scope':'260 ps pre-stretched charged-triplet no-weak '
                'long-time worldsheet stress test, 5 vertices only; '
                'not an annihilation or final NG prediction'}
    Path('results').mkdir(exist_ok=True)
    Path('results/worldsheet_longtime_stress.json').write_text(json.dumps(result,indent=2)+'\n')
    print('[longtime-worldsheet] '+json.dumps(result),flush=True)
    return result


if __name__=='__main__':
    run()
