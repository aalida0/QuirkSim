"""Rapidity + zero-interior-mass equal-arclength gauge, 2T long-time pilot.

The newer constrained semidiscrete Nambu-Goto string has no artificial
interior masses. This independently tests its velocity-coordinate domain
at large CM boosts by evolving w=gamma*v and retaining BOTH:
  - common total physical arclength stop,
  - separate local element-quality stops.
The result remains a gated numerical development test, not a physical
annihilation or full continuum worldsheet prediction.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

from quirksim.variational_string import StringParameters
from quirksim.equidistributed_ng import EquidistributedNGString
from quirksim.variational_rapidity import (
    physical_to_rapidity,rapidity_to_physical,rhs_rapidity)
from quirksim.mesh_guard import (
    make_total_arclength_stop_event,make_mesh_stop_events,
    total_string_arclength_cm)


def run(nodes=5,beta_cm_x=.65,horizon_ps=1000.,cutoff_fraction=1e-5):
    if not 0<cutoff_fraction<1:
        raise ValueError('invalid physical length cutoff fraction')
    p=StringParameters(nodes=nodes,initial_length_cm=.1,
                       mass_gev=225,lambda_ev=774,initial_beta=.4,
                       b_tesla=2.,gauge_mass_fraction=1e-3,
                       regularizer_cm=1e-9,max_step_ps=.5)
    engine=EquidistributedNGString(p)
    start_state=engine.boosted_straight_initial_state(beta_cm_x)
    initial=engine.action_quantities(start_state)
    w0=physical_to_rapidity(start_state,nodes)
    Lcut=p.initial_length_cm*cutoff_fraction
    stop=[make_total_arclength_stop_event(nodes,Lcut)]
    stop+=make_mesh_stop_events(nodes,p.regularizer_cm)
    try:
        sol=solve_ivp(lambda t,w:rhs_rapidity(engine,t,w),
             (0,horizon_ps),w0,method='DOP853',
             atol=1e-10,rtol=2e-9,max_step=p.max_step_ps,events=stop)
        if not sol.success:
            raise RuntimeError(sol.message)
        t=float(sol.t[-1])
        y=rapidity_to_physical(sol.y[:,-1],nodes)
        x,v=engine._unpack(y)
        final=engine.action_quantities(y)
        _,_,g,gdot=engine.constraints(y)
        triggers=[label for label,hits in zip(('total_arclength','short_element','element_ratio'),
                                               sol.t_events) if len(hits)]
        material=np.linalg.norm(sol.y[3*nodes:6*nodes,:].T.reshape(-1,nodes,3),axis=2)
        maxspeed=float(np.max(material/np.sqrt(1+material**2)))
        positions=sol.y[:3*nodes,:].T.reshape(-1,nodes,3)
        seg=np.linalg.norm(np.diff(positions,axis=1),axis=2)
        shortest=float(np.min(seg))
        worst=float(np.max(np.max(seg,axis=1)/
                           np.maximum(np.min(seg,axis=1),p.regularizer_cm)))
        return {
            'status':('LONGTIME_NUMERICAL_REACHED_HORIZON' if not triggers else
                      'STOPPED_AT_PHYSICAL_STRING_LENGTH' if triggers==['total_arclength'] else
                      'BLOCKED_MESH_QUALITY'),
            'triggered':triggers,'nodes':nodes,'beta_cm_x':beta_cm_x,
            'horizon_requested_ps':horizon_ps,'time_reached_ps':t,
            'physical_arclength_cutoff_cm':Lcut,
            'total_arclength_cm':total_string_arclength_cm(y,nodes),
            'endpoint_separation_cm':float(np.linalg.norm(x[0]-x[-1])),
            'max_vertex_beta':maxspeed,
            'min_segment_over_trajectory_cm':shortest,
            'worst_element_ratio_over_trajectory':worst,
            'gauge_length_violation_cm':float(np.max(np.abs(g))),
            'gauge_length_rate_violation_cm_per_ps':float(np.max(np.abs(gdot))),
            'energy_drift_gev':final['energy_gev']-initial['energy_gev'],
            'momentum_change_gev':float(np.linalg.norm(
                final['momentum_gev']-initial['momentum_gev'])),
            'nfev':sol.nfev,'gauge_mass_gev':engine.gauge_mass,
            'limits':['no physical radiation, weak decay, or annihilation',
                      'gated at string length cutoff; cannot infer collision probability',
                      'not a validated continuum worldsheet or exact lab boost']
        }
    except Exception as exc:
        return {'status':'BLOCKED_RAPIDITY_CONSTRAINED_NG',
                'nodes':nodes,'beta_cm_x':beta_cm_x,
                'horizon_requested_ps':horizon_ps,
                'reason':str(exc)}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--nodes',type=int,default=5)
    p.add_argument('--beta-cm-x',type=float,default=.65)
    p.add_argument('--horizon-ps',type=float,default=1000.)
    p.add_argument('--cutoff-fraction',type=float,default=1e-5)
    p.add_argument('--output',default='results/constrained_ng_rapidity.json')
    args=p.parse_args()
    result=run(args.nodes,args.beta_cm_x,args.horizon_ps,args.cutoff_fraction)
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+'\n')
    print('[constrained-rapidity] '+json.dumps(result),flush=True)


if __name__=='__main__':
    main()
