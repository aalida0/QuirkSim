"""Diagnose zero-interior-mass gauge-constrained Nambu--Goto first return.

Compares an equal-spatial-arclength gauge to the old artificial-mass
worldsheet near its first collapse. All states must pass mesh and Noether
gates; hitting the physical zero-length pinch is a distinct *blocked*
boundary problem, not a license to continue through it.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

from quirksim.equidistributed_ng import EquidistributedNGString
from quirksim.variational_string import StringParameters
from quirksim.mesh_guard import make_mesh_stop_events
from quirksim.radial_ng_reference import ExactRadialNG,RadialNGParameters


def case(n,b_tesla,beta_boost_x,maximum_ps=260.):
    p=StringParameters(nodes=n,mass_gev=225.,lambda_ev=774.,
                       initial_length_cm=.1,initial_beta=.4,
                       b_tesla=b_tesla,gauge_mass_fraction=1e-3,
                       regularizer_cm=1e-9,max_step_ps=.5)
    s=EquidistributedNGString(p)
    y=(s.initial_state() if beta_boost_x==0
       else s.boosted_straight_initial_state(beta_cm_x=beta_boost_x))
    start=s.action_quantities(y)
    ref=ExactRadialNG(RadialNGParameters(initial_separation_cm=.1,
                  mass_neutral_gev=225.,lambda_ev=774.,initial_beta=.4))
    try:
        sol=solve_ivp(s.rhs,(0,maximum_ps),y,method='DOP853',
            rtol=2e-9,atol=1e-10,max_step=p.max_step_ps,
            events=make_mesh_stop_events(n,p.regularizer_cm))
        if not sol.success:
            raise RuntimeError(sol.message)
        t=float(sol.t[-1])
        final=s.action_quantities(sol.y[:,-1])
        x,v=s._unpack(sol.y[:,-1])
        q=s.mesh_quality(sol.y)
        _,_,g,gdot=s.constraints(sol.y[:,-1])
        stops=[name for name, hits in zip(('short_element','element_ratio'),
                                      sol.t_events) if len(hits)]
        separation=float(np.linalg.norm(x[-1]-x[0]))
        comparable=b_tesla==0 and beta_boost_x==0
        reference_sep=ref.state(t)['unsigned_separation_cm']
        return {
          'status':'FIRST_MESH_STOP' if stops else 'REACHED_HORIZON',
          'stops':stops,'nodes':n,'b_tesla':b_tesla,
          'beta_boost_x':beta_boost_x,'time_ps':t,
          'time_divided_by_radial_first_crossing':t/ref.first_crossing_ps,
          'endpoint_separation_cm':separation,
          'exact_radial_reference_separation_cm':reference_sep,
          'exact_radial_reference_comparable':comparable,
          'radial_reference_absolute_separation_error_cm':
              abs(separation-reference_sep) if comparable else None,
          'max_equal_length_constraint_cm':float(np.max(np.abs(g))),
          'max_equal_length_rate_constraint_cm_per_ps':
              float(np.max(np.abs(gdot))),
          'mesh':q,
          'gauge_point_mass_gev':s.gauge_mass,
          'energy_drift_gev':final['energy_gev']-start['energy_gev'],
          'momentum_drift_gev':float(np.linalg.norm(
             final['momentum_gev']-start['momentum_gev'])),
          'rhs_evaluations':sol.nfev,
          'physics_status':'NOT_CONTINUUM_VALIDATED_NO_POST_PINCH_EVOLUTION',
        }
    except Exception as exc:
        return {'status':'BLOCKED_BY_CONSTRAINED_GAUGE_SOLVER',
                'nodes':n,'b_tesla':b_tesla,
                'beta_boost_x':beta_boost_x,'reason':str(exc)}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--nodes',type=int,default=5)
    p.add_argument('--max-time-ps',type=float,default=260.)
    p.add_argument('--output',default='results/equidistributed_gate.json')
    args=p.parse_args()
    records=[]
    for b,boost in ((0.,0.),(0.,.25),(2.,.25)):
        result=case(args.nodes,b,boost,args.max_time_ps)
        records.append(result)
        print('[equidistributed-gauge] '+json.dumps(result),flush=True)
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps({'records':records,
        'scope':'zero-mass-interior gauge-constrained finite-element test, '
                'no weak decays or physical crossing continuation'},indent=2)+'\n')


if __name__=='__main__':
    main()
