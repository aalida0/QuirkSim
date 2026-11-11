"""Rapidity-coordinate integrator for variational string, |v|<c by design.

For every vertex w=gamma*v, so v=w/sqrt(1+w.w). The differential
transformation is dw/dt=gamma*(dv/dt+gamma^2*v*(v.dv/dt)).
The semidiscrete Nambu-Goto ACTION and forces are not modified.
This may prevent accidental superluminal *intermediate RK stages*.
It does NOT cure poorly resolved folds, gauge modes, or continuum issues.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from .mesh_guard import make_mesh_stop_events,require_good_mesh,assert_no_mesh_stop


def physical_to_rapidity(y,nodes):
    out=np.array(y,float,copy=True)
    v=out[3*nodes:6*nodes].reshape(nodes,3)
    v2=np.einsum('ij,ij->i',v,v)
    if np.max(v2)>=1:
        raise ValueError('physical_to_rapidity requires subluminal input')
    v*=1/np.sqrt(1-v2)[:,None]
    return out


def rapidity_to_physical(y,nodes):
    out=np.array(y,float,copy=True)
    w=out[3*nodes:6*nodes].reshape(nodes,3)
    w2=np.einsum('ij,ij->i',w,w)
    w*=1/np.sqrt(1+w2)[:,None]
    return out


def rhs_rapidity(engine,t,yw):
    n=engine.n
    yp=rapidity_to_physical(yw,n)
    dp=engine.rhs(t,yp)
    w=yw[3*n:6*n].reshape(n,3)
    w2=np.einsum('ij,ij->i',w,w)
    gamma=np.sqrt(1+w2)
    v=yp[3*n:6*n].reshape(n,3)
    dv=dp[3*n:6*n].reshape(n,3)
    dw=gamma[:,None]*(dv+gamma[:,None]**2*
                       v*np.sum(v*dv,axis=1)[:,None])
    out=dp.copy()
    out[3*n:6*n]=dw.ravel()
    return out


def solve_rapidity(engine,state,duration_ps,max_step_ps=None,strict_mesh=True):
    n=engine.n
    if max_step_ps is None:
        max_step_ps=engine.params.max_step_ps
    if strict_mesh:
        require_good_mesh(state,n,engine.params.regularizer_cm)
    start=engine.action_quantities(state)
    w0=physical_to_rapidity(state,n)
    mesh_events=(make_mesh_stop_events(n,engine.params.regularizer_cm)
                 if strict_mesh else None)
    sol=solve_ivp(lambda t,y:rhs_rapidity(engine,t,y),
        (0,duration_ps),w0,method='DOP853',atol=1e-10,rtol=2e-9,
        max_step=max_step_ps,events=mesh_events)
    if not sol.success:
        raise RuntimeError(f'rapidity-string solver: {sol.message}')
    if strict_mesh:
        assert_no_mesh_stop(sol)
    velocities=sol.y[3*n:6*n,:].T.reshape(-1,n,3)
    vmax=float(np.max(np.sqrt(np.sum(velocities**2,axis=2)/
                            (1+np.sum(velocities**2,axis=2)))))
    max_w=float(np.max(np.linalg.norm(velocities,axis=2)))
    sol.y=np.stack([rapidity_to_physical(row,n) for row in sol.y.T],axis=1)
    end=engine.action_quantities(sol.y[:,-1])
    return {
       'initial_energy_gev':start['energy_gev'],
       'final_energy_gev':end['energy_gev'],
       'energy_drift_gev':end['energy_gev']-start['energy_gev'],
       'momentum_drift_gev':float(np.linalg.norm(
           end['momentum_gev']-start['momentum_gev'])),
       'angular_momentum_drift_over_hbar':float(np.linalg.norm(
           end['angular_momentum_over_hbar']-
           start['angular_momentum_over_hbar'])),
       'nfev':sol.nfev,
       'max_vertex_beta':vmax,
       'max_vertex_spatial_rapidity':max_w,
       'caveat':'subluminal-coordinate solver does not prove continuum causality'
    },sol
