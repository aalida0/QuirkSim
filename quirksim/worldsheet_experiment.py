"""Proper-time clock pilot on the variational NG finite-element string.

Both endpoints begin electrically charged and change to neutral after
Q+/- -> Q0 + pion / charged lepton + neutrino. The finite-element weak
impulse is a conservative canonical momentum update with a STRICT energy
gate. If a physical gate or ODE fails, return BLOCKED, not a fake event.

Research limits: initial pre-stretched lab-time slice rather than pointlike
production, gauge mass regulator and spatial quadrature not continuum
converged, no radiation/annihilation or string self-intersection model.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

from .dynamics import C_CM_PER_PS,sample_plan
from .physics import benchmark
from .variational_string import StringParameters,VariationalString
from .encounter_hazards import scenario_envelope


def pair_boost_initial_state(s,beta_cm_x):
    """Simplified lab slice of a pre-stretched string; not an exact Lorentz
    transformation of a history-dependent moving Nambu-Goto worldsheet."""
    if abs(beta_cm_x)>=1:
        raise ValueError('invalid pair boost')
    y=s.initial_state()
    x,v=s._unpack(y)
    beta_rel=s.params.initial_beta
    v[0,0]=(-beta_rel+beta_cm_x)/(1-beta_rel*beta_cm_x)
    v[-1,0]=(beta_rel+beta_cm_x)/(1+beta_rel*beta_cm_x)
    v[1:-1,0]=beta_cm_x
    x[:,0]*=math.sqrt(1-beta_cm_x**2)
    return y


def simulate(seed=25,params=None,duration_ps=5, beta_cm_x=.25,
             energy_tolerance_gev=1e-4):
    if params is None:
        params=StringParameters(nodes=5,initial_length_cm=.1,
              initial_beta=.4,b_tesla=2,max_step_ps=.05)
    engine=VariationalString(params)
    bench=benchmark('triplet',params.mass_gev)
    latent=sample_plan(bench['lifetime_ps'],bench['widths_gev'],seed)
    y=pair_boost_initial_state(engine,beta_cm_x)
    t=0.
    events=[]
    closest_minima=[]
    segment_drift=[]
    failed=None
    while t<duration_ps-1e-11:
        active=[i for i in (0,1) if engine.charge[[0,engine.n-1][i]]!=0]
        if not active:
            active=[]
        n=engine.n

        def radial_min(_t,state):
            x,v=engine._unpack(state)
            r=x[0]-x[-1]
            return float(r@(v[0]-v[-1]))
        radial_min.direction=1
        radial_min.terminal=False

        def decay_event(index):
            def e(_t,state):
                return state[-2+index]-latent[index]['proper_time_ps']
            e.direction=1
            e.terminal=True
            return e

        initial_E=engine.action_quantities(y)['energy_gev']
        try:
            sol=solve_ivp(engine.rhs,(t,duration_ps),y,method='DOP853',
                          atol=1e-10,rtol=2e-9,max_step=params.max_step_ps,
                          events=[radial_min]+[decay_event(i) for i in active])
        except (RuntimeError,ValueError,np.linalg.LinAlgError) as exc:
            failed=f'worldsheet integration failed: {exc}'
            break
        if not sol.success:
            failed=f'worldsheet integrator: {sol.message}'
            break
        for tt,zz in zip(sol.t_events[0],sol.y_events[0]):
            if tt>1e-6:
                x,_=engine._unpack(zz)
                closest_minima.append({'time_ps':float(tt),
                    'separation_cm':float(np.linalg.norm(x[0]-x[-1]))})
        t=float(sol.t[-1])
        y=sol.y[:,-1].copy()
        segment_drift.append(float(engine.action_quantities(y)['energy_gev']-initial_E))
        triggered=[i for k,i in enumerate(active) if len(sol.t_events[k+1])]
        if not triggered:
            break
        label=triggered[0]
        endpoint=[0,engine.n-1][label]
        channel=latent[label]['channel']
        try:
            y,info=engine.weak_decay_impulse(
                y,endpoint,channel,seed=latent[label]['rest_seed'],
                energy_tolerance_gev=energy_tolerance_gev,strict=True)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            failed=f'weak-decay canonical energy gate: {exc}'
            break
        events.append({'time_ps':t,'endpoint_label':label,
                       'proper_time_ps':float(y[-2+label]),**info})
    first=events[0]['time_ps'] if events else None
    post=[e for e in closest_minima if first is not None and
          e['time_ps']>first+1e-7]
    first_hazard=scenario_envelope(closest_minima,first,
                                   params.lambda_ev,params.mass_gev)
    r={
      'status':'BLOCKED_PHYSICS_GATE' if failed else 'FINITE_ELEMENT_PROTOTYPE_ONLY',
      'failure_reason':failed,'seed':seed,'parameters':vars(params),
      'pair_boost_x':beta_cm_x,'elapsed_ps':t,
      'charged_endpoints':[bool(engine.charge[0]),bool(engine.charge[-1])],
      'latent_decay_plan':latent,
      'decays':events,
      'string_energy_drift_by_free_segment_gev':segment_drift,
      'geometric_closest_encounters':closest_minima,
      'post_decay_closest_cm':min((e['separation_cm'] for e in post),default=None),
      'radiation_annihilation_hazard_envelope_summary':{
        'scenarios':len(first_hazard),
        'possible_survival_to_weak_range':[
           min((x['conditional_survival_to_first_weak_decay'] for x in first_hazard
                if x['conditional_survival_to_first_weak_decay'] is not None),default=None),
           max((x['conditional_survival_to_first_weak_decay'] for x in first_hazard
                if x['conditional_survival_to_first_weak_decay'] is not None),default=None)],
        'note':'fully assumed contact hazard scenarios; NO self-consistent emission'},
      'physics_limits':[
         'not exact continuum NG worldsheet; mesh and gauge regulator',
         'stretched-string startup and approximate lab-frame boost',
         'no physical annihilation survival or radiation recoil in evolution',
         'isolated piloted weak events do not define inclusive pp predictions']
      }
    return r


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=int,default=25)
    p.add_argument('--duration-ps',type=float,default=5)
    p.add_argument('--nodes',type=int,default=5)
    p.add_argument('--b-tesla',type=float,default=2)
    p.add_argument('--beta-cm-x',type=float,default=.25)
    p.add_argument('--gauge-mass-fraction',type=float,default=.001)
    p.add_argument('--energy-tolerance-gev',type=float,default=1e-4)
    p.add_argument('--output',default='results/variational_weak_pilot.json')
    a=p.parse_args()
    params=StringParameters(nodes=a.nodes,b_tesla=a.b_tesla,
          gauge_mass_fraction=a.gauge_mass_fraction,
          max_step_ps=.05)
    result=simulate(seed=a.seed,params=params,duration_ps=a.duration_ps,
           beta_cm_x=a.beta_cm_x,energy_tolerance_gev=a.energy_tolerance_gev)
    file=Path(a.output);file.parent.mkdir(parents=True,exist_ok=True)
    file.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
