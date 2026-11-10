"""Convergence tests of magnetic, recoil-driven nonradial quirk trajectories.

Numerical stability is not physics-model validation. This checks the adaptive
ODE step and short-distance regularization for nonzero post-decay misses.
"""
import json
from dataclasses import replace
from pathlib import Path

from quirksim.field_dynamics import FieldParameters, simulate


def run():
    outcomes=[]
    for seed in (10,42,197):
        base=FieldParameters(seed=seed,model='triplet',mass_gev=225,
                             lambda_ev=774,b_tesla=2,beta_cm_x=.25,
                             beta_rel=.4,mode='recoil',force_model='ng_endpoint',
                             duration_ps=1000,max_step_ps=2.,regulator_cm=1e-7)
        a=simulate(base)
        b=simulate(replace(base,max_step_ps=.5))
        c=simulate(replace(base,regulator_cm=5e-8,max_step_ps=.5))
        assert a['plan']==b['plan']==c['plan']
        assert [e['channel'] for e in a['decays']]==[e['channel'] for e in b['decays']]==[e['channel'] for e in c['decays']]
        d=[x['closest_post_decay_turning_point_cm'] for x in (a,b,c)]
        if any(x is None for x in d):
            raise RuntimeError('undefined closest post-decay encounter')
        dr_step=1e4*abs(d[0]-d[1])
        dr_reg=1e4*abs(d[1]-d[2])
        row={'seed':seed,'rmin_um_default':d[0]*1e4,
             'rmin_um_small_step':d[1]*1e4,
             'rmin_um_half_regulator':d[2]*1e4,
             'abs_solver_step_diff_um':dr_step,
             'abs_regulator_diff_um':dr_reg}
        print('[convergence] '+json.dumps(row),flush=True)
        outcomes.append(row)
    summary={'rows':outcomes,
             'max_solver_step_diff_um':max(x['abs_solver_step_diff_um'] for x in outcomes),
             'max_regulator_diff_um':max(x['abs_regulator_diff_um'] for x in outcomes),
             'scope':'3 selected recoil seeds, NG straight-string endpoint-force approximation, 2T solenoid; not full worldsheet validation'}
    Path('results').mkdir(exist_ok=True)
    Path('results/triplet_field_convergence.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('[convergence-summary] '+json.dumps(summary),flush=True)
    if summary['max_solver_step_diff_um']>.05 or summary['max_regulator_diff_um']>.05:
        raise AssertionError('nonradial numerical convergence check failed')
    return summary


if __name__=='__main__':
    run()
