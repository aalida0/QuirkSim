"""Bounded-speed rapidity coordinates + zero-mass constrained NG spatial gauge.

This combines the constrained-gauge action solver with w=gamma*v variables
(there is no change to the physical semidiscrete NG Lagrangian).  No 3D
contact continuation, string reconnection or annihilation is implemented.

An event stops the evolution at the first crossing of a total-arclength
cutoff, while independent guards stop invalid finite elements.  The same
geometric cutoff is used for every spatial mesh resolution.
"""
from dataclasses import dataclass
import json
from pathlib import Path
import math
import numpy as np
from scipy.integrate import solve_ivp

from quirksim.equidistributed_ng import EquidistributedNGString
from quirksim.variational_string import StringParameters
from quirksim.variational_rapidity import (
    physical_to_rapidity,rapidity_to_physical,rhs_rapidity)
from quirksim.mesh_guard import (
    make_mesh_stop_events,make_total_arclength_stop_event,
    total_string_arclength_cm)
from quirksim.dynamics import HBARC_GEV_CM


def run_case(nodes=11, b_tesla=2., beta_cm_x=.25,
             duration_ps=1000., pinch_fraction=1e-6,
             regulator_cm=1e-12, max_step_ps=.2,
             rtol=2e-9,atol=1e-10,method='DOP853'):
    params=StringParameters(nodes=nodes,mass_gev=225.,lambda_ev=774.,
        initial_length_cm=.1,initial_beta=.4,
        b_tesla=b_tesla,regularizer_cm=regulator_cm,
        gauge_mass_fraction=1e-3,max_step_ps=max_step_ps)
    s=EquidistributedNGString(params)
    y=s.boosted_straight_initial_state(beta_cm_x=beta_cm_x)
    pre=s.action_quantities(y)
    w0=physical_to_rapidity(y,nodes)
    cutoff=pinch_fraction*params.initial_length_cm
    events=[make_total_arclength_stop_event(nodes,cutoff)]
    events+=make_mesh_stop_events(nodes,regulator_cm)
    try:
        if method not in ('DOP853','Radau','BDF','LSODA'):
            raise ValueError('unsupported ODE method')
        sol=solve_ivp(lambda t,w:rhs_rapidity(s,t,w),(0,duration_ps),w0,
            method=method,rtol=rtol,atol=atol,
            max_step=max_step_ps,events=events)
        if not sol.success:
            last=rapidity_to_physical(sol.y[:,-1],nodes)
            lastx,lastv=s._unpack(last)
            return {
                'status':'BLOCKED_ODE_STIFFNESS',
                'nodes':nodes,'method':method,'reason':sol.message,
                'last_valid_time_ps':float(sol.t[-1]),
                'last_valid_total_string_arclength_cm':
                    total_string_arclength_cm(last,nodes),
                'last_valid_endpoint_chord_cm':
                    float(np.linalg.norm(lastx[-1]-lastx[0])),
                'last_valid_max_vertex_beta':
                    float(np.max(np.linalg.norm(lastv,axis=1))),
                'infracolor_reference_length_cm':HBARC_GEV_CM/(774e-9),
                'requested_arclength_cutoff_cm':cutoff,
                'physics_gate':'BLOCKED_NO_POST_FAILURE_CONTINUATION'
            }
        names=('common_arclength_cutoff','short_element','element_imbalance')
        stops=[name for name,e in zip(names,sol.t_events) if len(e)]
        final=rapidity_to_physical(sol.y[:,-1],nodes)
        after=s.action_quantities(final)
        positions,velocity=s._unpack(final)
        gauge=s.constraints(final)
        length=total_string_arclength_cm(final,nodes)
        chord=float(np.linalg.norm(positions[-1]-positions[0]))
        dv=np.diff(positions,axis=0)
        tangent=dv/np.linalg.norm(dv,axis=1)[:,None]
        adjacent=np.arccos(np.clip(np.sum(tangent[:-1]*tangent[1:],axis=1),-1.,1.))
        maximum_beta=float(np.max(np.linalg.norm(velocity,axis=1)))
        gauge_residual=float(np.max(np.abs(gauge[2])))
        reference=HBARC_GEV_CM/(774e-9)
        return {
            'status':('REACHED_COMMON_TOTAL_ARCLENGTH' if
                      stops==['common_arclength_cutoff'] else
                      'REACHED_HORIZON' if not stops else
                      'BLOCKED_PREMATURE_MESH_STOP'),
            'stops':stops,
            'nodes':nodes,'b_tesla':b_tesla,'beta_cm_x':beta_cm_x,
            'duration_requested_ps':duration_ps,
            'time_ps':float(sol.t[-1]),
            'regulator_cm':regulator_cm,'rtol':rtol,'atol':atol,
            'method':method,
            'max_step_ps':max_step_ps,
            'physical_total_arclength_cm':length,
            'common_total_arclength_cutoff_cm':cutoff,
            'arclength_event_residual_cm':abs(length-cutoff)
                if 'common_arclength_cutoff' in stops else None,
            'endpoint_chord_cm':chord,
            'chord_deficit_cm':length-chord,
            'maximum_adjacent_bend_angle_rad':float(np.max(adjacent)),
            'max_vertex_beta':maximum_beta,
            'min_segment_length_cm':float(np.min(np.linalg.norm(dv,axis=1))),
            'max_gauge_length_deviation_cm':gauge_residual,
            'energy_drift_gev':after['energy_gev']-pre['energy_gev'],
            'mechanical_momentum_change_gev':float(np.linalg.norm(
                after['momentum_gev']-pre['momentum_gev'])),
            'chord_deficit_in_1_over_lambda':(length-chord)/reference,
            'nfev':sol.nfev,
            'physics_scope':'bounded-speed coordinate choice, zero interior masses; '
                 'NOT validated continuum NG, post-contact, radiation or annihilation'
        }
    except Exception as exc:
        return {'status':'BLOCKED_RAPIDITY_CONSTRAINED_SOLVER','nodes':nodes,
                'regulator_cm':regulator_cm,'requested_cutoff_cm':cutoff,
                'reason':str(exc)}


def main():
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('--nodes',type=int,default=11)
    p.add_argument('--regulator-cm',type=float,default=1e-12)
    p.add_argument('--pinch-fraction',type=float,default=1e-6)
    p.add_argument('--max-step-ps',type=float,default=.2)
    p.add_argument('--method',choices=['DOP853','Radau','BDF','LSODA'],default='DOP853')
    p.add_argument('--duration-ps',type=float,default=1000.)
    p.add_argument('--output',default='results/rapidity_constrained_contact.json')
    args=p.parse_args()
    if not 0<args.pinch_fraction<1:
        p.error('--pinch-fraction must be (0,1)')
    r=run_case(nodes=args.nodes,regulator_cm=args.regulator_cm,
         pinch_fraction=args.pinch_fraction,max_step_ps=args.max_step_ps,
         duration_ps=args.duration_ps,method=args.method)
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(r,indent=2)+'\n')
    print('[rapidity-constrained-contact] '+json.dumps(r),flush=True)


if __name__=='__main__':
    main()
