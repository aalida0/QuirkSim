"""Compare transverse string vibrations to the linearized NG wave equation.

For a stretched static string, small transverse fluctuations obey
d²y/dt² = c² d²y/dx². The lowest fixed-end mode has a half-period L/c.
Our heavy but finite endpoints slightly change the continuum reference;
this is a refinement/stability *diagnostic*, not a proof of convergence.
"""
import json
from pathlib import Path
import math
import numpy as np
from quirksim.variational_string import StringParameters,VariationalString
from quirksim.dynamics import C_CM_PER_PS


def run():
    length=.05
    amp=1e-5
    half_period_ps=length/C_CM_PER_PS
    records=[]
    for n in (5,7,11):
        p=StringParameters(nodes=n,initial_beta=0.,
                            initial_length_cm=length,duration_ps=half_period_ps,
                            gauge_mass_fraction=1e-3,max_step_ps=.04)
        s=VariationalString(p)
        state=s.initial_state()
        x,v=s._unpack(state)
        x[:,1]=amp*np.sin(np.pi*np.linspace(0,1,n))
        e_before=s.action_quantities(state)['energy_gev']
        diag,sol=s.solve(state)
        qfinal,vfinal=s._unpack(sol.y[:,-1])
        center_ratio=float(qfinal[n//2,1]/amp)
        row={'nodes':n,'half_period_ps':half_period_ps,
             'center_y_over_initial_after_half_period':center_ratio,
             'distance_from_linear_NG_prediction_abs':abs(center_ratio+1.),
             'total_energy_drift_gev':diag['energy_drift_gev'],
             'nfev':diag['n_rhs_evaluations']}
        print('[wave-mode] '+json.dumps(row),flush=True)
        records.append(row)
    report={'records':records,
            'notes':'Lowest transverse linear mode with finite massive endpoints and '
                    'small tangent gauge mass. Decreasing error under spatial refinement '
                    'is necessary but not sufficient for full-worldsheet validation.',
            'refinement_improves_half_period_error':
                records[-1]['distance_from_linear_NG_prediction_abs']<
                records[0]['distance_from_linear_NG_prediction_abs']}
    Path('results').mkdir(exist_ok=True)
    Path('results/variational_transverse_wave.json').write_text(json.dumps(report,indent=2)+'\n')
    print('[wave-summary] '+json.dumps(report),flush=True)
    return report


if __name__=='__main__':
    run()
