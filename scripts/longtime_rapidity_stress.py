"""Rapidity-worldsheet 260 ps stress with explicit mesh-quality diagnostics."""
import argparse
import json
import numpy as np
from pathlib import Path
from quirksim.variational_string import StringParameters,VariationalString
from quirksim.variational_rapidity import solve_rapidity
from quirksim.worldsheet_experiment import pair_boost_initial_state


def run(b_field=0,duration=260,nodes=5):
    p=StringParameters(nodes=nodes,mass_gev=225,lambda_ev=774,
                       initial_length_cm=.1,initial_beta=.4,
                       b_tesla=b_field,max_step_ps=.5,
                       gauge_mass_fraction=1e-3,duration_ps=duration)
    s=VariationalString(p)
    y=pair_boost_initial_state(s,.25)
    try:
        diag,sol=solve_rapidity(s,y,duration,.5)
        pos=sol.y[:3*nodes,:].T.reshape(-1,nodes,3)
        edges=np.linalg.norm(np.diff(pos,axis=1),axis=2)
        short=float(edges.min())
        worst=float(np.max(edges.max(axis=1)/np.maximum(edges.min(axis=1),p.regularizer_cm)))
        start=s.action_quantities(y);end=s.action_quantities(sol.y[:,-1])
        report={'status':'NUMERICAL_RUN_COMPLETE_MESH_GATED',
            'b_field_tesla':b_field,'n':nodes,'duration_ps':duration,
            'energy_drift_gev':diag['energy_drift_gev'],
            'momentum_drift_gev':diag['momentum_drift_gev'],
            'max_vertex_beta':diag['max_vertex_beta'],
            'max_vertex_rapidity':diag['max_vertex_spatial_rapidity'],
            'min_string_element_cm':short,'worst_length_ratio':worst,
            'mesh_collapsed':short<100*p.regularizer_cm,
            'mesh_imbalanced':worst>50,
            'nfev':sol.nfev,
            'physics_gate':'BLOCKED' if short<100*p.regularizer_cm or worst>50
                 or abs(diag['energy_drift_gev'])>1e-3 else 'LONGTIME_NUMERICAL_PASS_ONLY',
            'note':'even a numerical pass is NOT a continuum worldsheet solution'}
    except Exception as exc:
        report={'status':'BLOCKED_BY_NUMERICAL_STRING_EVOLUTION',
                'b_field_tesla':b_field,'n':nodes,
                'duration_ps_requested':duration,'reason':str(exc)}
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--b-tesla',type=float,default=0.)
    p.add_argument('--duration-ps',type=float,default=260.)
    p.add_argument('--nodes',type=int,default=5)
    p.add_argument('--output',default='results/longtime_rapidity.json')
    args=p.parse_args()
    result=run(args.b_tesla,args.duration_ps,args.nodes)
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n')
    print('[rapidity-stress] '+json.dumps(result),flush=True)


if __name__=='__main__':
    main()
