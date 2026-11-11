"""Pre-pinch convergence of constrained zero-gauge-mass Nambu--Goto string.

Runs the same pre-stretched triplet B=2T boosted trajectory with 5,7,11
vertices to a COMMON physical time 200 ps (before the first close return).
No weak decays, string radiation or annihilation. Compare end positions
and energies, not gauge vertex trajectories.
"""
import json
from pathlib import Path
import numpy as np

from quirksim.variational_string import StringParameters
from quirksim.equidistributed_ng import EquidistributedNGString


def run():
    rows=[]
    for n in (5,7,11):
        p=StringParameters(nodes=n,mass_gev=225,lambda_ev=774,
                           initial_length_cm=.1,initial_beta=.4,
                           b_tesla=2.,regularizer_cm=1e-9,
                           max_step_ps=.5)
        s=EquidistributedNGString(p)
        state=s.boosted_straight_initial_state(beta_cm_x=.25)
        d,sol=s.solve(state,duration_ps=200,
                      pinch_fraction_of_initial_length=4e-4)
        x,v=s._unpack(sol.y[:,-1])
        row={
            'n':n,'time_ps':200.,
            'endpoint_left_cm':x[0].tolist(),
            'endpoint_right_cm':x[-1].tolist(),
            'endpoint_separation_cm':float(np.linalg.norm(x[0]-x[-1])),
            'endpoint_left_velocity_c':v[0].tolist(),
            'endpoint_right_velocity_c':v[-1].tolist(),
            'energy_drift_gev':float(d['energy_drift_gev']),
            'max_equal_arclength_deviation_cm':
                d['max_equal_arclength_deviation_cm'],
            'mesh_quality':d['mesh_quality']
        }
        rows.append(row)
        print('[prepinch] '+json.dumps(row),flush=True)
    final=rows[-1]
    for row in rows[:-1]:
        row['position_difference_from_11node_max_cm']=max(
            float(np.linalg.norm(np.array(row['endpoint_left_cm'])-
                  np.array(final['endpoint_left_cm']))),
            float(np.linalg.norm(np.array(row['endpoint_right_cm'])-
                  np.array(final['endpoint_right_cm']))))
    output={'status':'PREPINCH_MESH_REFINEMENT_ONLY',
            'rows':rows,
            'largest_5_7_position_difference_vs_11_cm':max(
                row['position_difference_from_11node_max_cm'] for row in rows[:-1]),
            'limits':'200ps only; no worldsheet folding, weak recoil, '
                     'string emission or collision continuation'}
    Path('results').mkdir(exist_ok=True)
    Path('results/equidistributed_prepinch_convergence.json').write_text(
        json.dumps(output,indent=2)+'\n')
    print('[prepinch-summary] '+json.dumps({
        'largest_position_diff_cm':output['largest_5_7_position_difference_vs_11_cm'],
        'max_energy_drift_gev':max(abs(r['energy_drift_gev']) for r in rows)}),
        flush=True)
    if (output['largest_5_7_position_difference_vs_11_cm']>1e-4
            or max(abs(r['energy_drift_gev']) for r in rows)>1e-4):
        raise AssertionError('constrained NG pre-pinch mesh refinement gate')
    return output


if __name__=='__main__':
    run()
