"""First-violation diagnostic for finite-element string mesh/rapidity solver.

Stop at FIRST underresolved element or 50:1 element imbalance; do not
integrate across a collapsed material coordinate. Compare B=0, unboosted
radial trajectory to the independent exact 1D massive-endpoint benchmark.
This is a diagnosis; it is NOT a remeshing or physics restoration.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

from quirksim.variational_string import StringParameters,VariationalString
from quirksim.variational_rapidity import (physical_to_rapidity,
    rapidity_to_physical,rhs_rapidity)
from quirksim.worldsheet_experiment import pair_boost_initial_state
from quirksim.radial_ng_reference import ExactRadialNG,RadialNGParameters


def diagnose(b_tesla=0.,beta_boost_x=0.,nodes=5,
             gauge_mass_fraction=1e-3,maximum_ps=260.,max_step_ps=.5):
    p=StringParameters(nodes=nodes,mass_gev=225.,lambda_ev=774.,
        initial_length_cm=.1,initial_beta=.4,
        b_tesla=b_tesla,gauge_mass_fraction=gauge_mass_fraction,
        max_step_ps=max_step_ps)
    engine=VariationalString(p)
    initial=(engine.initial_state() if beta_boost_x==0
             else pair_boost_initial_state(engine,beta_boost_x))
    w=physical_to_rapidity(initial,nodes)
    initial_conserved=engine.action_quantities(initial)
    reference=ExactRadialNG(RadialNGParameters(
        mass_neutral_gev=p.mass_gev,lambda_ev=p.lambda_ev,
        initial_beta=p.initial_beta,
        initial_separation_cm=p.initial_length_cm))

    def lengths(y):
        pos=np.asarray(y[:3*nodes]).reshape(nodes,3)
        return np.linalg.norm(np.diff(pos,axis=0),axis=1)

    def min_length_event(t,y):
        return float(np.min(lengths(y)) - 100*p.regularizer_cm)
    min_length_event.terminal=True
    min_length_event.direction=-1

    def imbalance_event(t,y):
        ell=lengths(y)
        return float(50. - ell.max()/max(ell.min(),p.regularizer_cm))
    imbalance_event.terminal=True
    imbalance_event.direction=-1

    try:
        sol=solve_ivp(lambda t,y:rhs_rapidity(engine,t,y),
           (0,maximum_ps),w,method='DOP853',max_step=max_step_ps,
           rtol=2e-9,atol=1e-10,
           events=[min_length_event,imbalance_event])
        if not sol.success:
            raise RuntimeError(sol.message)
        t=float(sol.t[-1])
        final=rapidity_to_physical(sol.y[:,-1],nodes)
        q=engine.action_quantities(final)
        ell=lengths(sol.y[:,-1])
        xs=final[:3*nodes].reshape(nodes,3)
        maxbeta=float(np.max(np.linalg.norm(
            final[3*nodes:6*nodes].reshape(nodes,3),axis=1)))
        fired=['short_element','element_ratio']
        triggers=[fired[k] for k in range(2) if len(sol.t_events[k])]
        rel_dist=float(np.linalg.norm(xs[-1]-xs[0]))
        radial_at_t=reference.state(t)
        report={
            'status':'MESH_BREAKDOWN_DETECTED' if triggers else
                     'OBSERVATION_END_NO_MESH_TRIGGER',
            'triggered':triggers,'time_ps':t,
            'b_tesla':b_tesla,'beta_boost_x':beta_boost_x,
            'nodes':nodes,'gauge_mass_fraction':gauge_mass_fraction,
            'minimum_element_cm':float(ell.min()),
            'maximum_element_cm':float(ell.max()),
            'element_imbalance_ratio':float(ell.max()/ell.min()),
            'endpoint_separation_cm':rel_dist,
            'max_vertex_beta':maxbeta,
            'total_energy_drift_gev':
                float(q['energy_gev']-initial_conserved['energy_gev']),
            'total_linear_momentum_drift_gev':
                float(np.linalg.norm(q['momentum_gev']-
                          initial_conserved['momentum_gev'])),
            'radial_reference_first_crossing_ps':
                reference.first_crossing_ps,
            'time_over_radial_first_crossing':
                t/reference.first_crossing_ps,
            'exact_unboosted_B0_radial_separation_cm_at_t':
                radial_at_t['unsigned_separation_cm'],
            'radial_reference_is_comparable':
                bool(b_tesla==0 and beta_boost_x==0),
            'nfev':sol.nfev,
            'physical_scope':'first mesh failure; not proof of continuum string or radiation/annihilation'
        }
    except Exception as exc:
        report={'status':'SOLVER_FAILED_BEFORE_MESH_EVENT',
                'b_tesla':b_tesla,'beta_boost_x':beta_boost_x,
                'nodes':nodes,'gauge_mass_fraction':gauge_mass_fraction,
                'reason':str(exc)}
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='results/mesh_first_failure.json')
    parser.add_argument('--max-time-ps',type=float,default=260.)
    args=parser.parse_args()
    studies=[]
    for nodes in (5,7):
        for field,boost in ((0.,0.),(0.,.25),(2.,.25)):
            record=diagnose(b_tesla=field,beta_boost_x=boost,
                            nodes=nodes,maximum_ps=args.max_time_ps)
            studies.append(record)
            print('[first-mesh-violation] '+json.dumps(record),flush=True)
    obj={'records':studies,
         'meaning':'gauge/mesh diagnostic compared to exact reduced radial NG '
                   'reference ONLY for unboosted B=0 case',
         'fail_closed':True}
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2)+'\n')
    return obj


if __name__=='__main__':
    main()
