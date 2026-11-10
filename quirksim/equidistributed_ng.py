"""Experimental equal-arclength gauge for an action-derived NG string.

Unlike VariationalString, which assigns UNPHYSICAL point masses to internal
nodes, this formulation sets interior masses to ZERO and imposes N-2
equal-spatial-segment-length gauge constraints. Their time derivatives fix the
null tangential components of the semidiscrete mass matrix.

Constraint:
    g_j = |x[j+1]-x[j]| - |x[j]-x[j-1]| = 0, 1<=j<N-1
    gddot_j = J_j @ dv/dt + b_j = 0.
Solve the variational saddle system
    [ M J.T ] [ dv/dt ] = [ physical_force - convective ]
    [ J  0  ] [ lambda ]   [ -b                      ].

Gauge forces are orthogonal to allowed velocity variations (J v = 0).
If numerical gauge drift becomes large, fail closed: do not extrapolate.
The action is a MIDPOINT spatial discretization and massive-endpoint
worldsheet self-intersection, finite propagation, and continuum convergence
are not fully solved by this gauge choice.

This experimental solver deliberately never switches to the old positive
artificial gauge masses.
"""
from dataclasses import replace
import math
import numpy as np
from scipy.integrate import solve_ivp

from .dynamics import C_CM_PER_PS
from .variational_string import VariationalString
from .mesh_guard import (make_mesh_stop_events,MeshQualityViolation,
                         require_good_mesh,assert_no_mesh_stop)


class EquidistributedNGString(VariationalString):
    def __init__(self,params):
        super().__init__(params)
        self.masses[1:-1]=0.
        self.gauge_mass=0.
        self.gauge_name='equal_spatial_arclength_zeros_interior_point_masses'

    def initial_state(self,transverse_velocity=0.,middle_offset_cm=0.):
        if middle_offset_cm:
            raise ValueError('bent initial geometry needs gauge-compatible velocity projection')
        y=super().initial_state(transverse_velocity=transverse_velocity)
        x,v=self._unpack(y)
        for i in range(1,self.n-1):
            u=i/(self.n-1)
            v[i]=(1-u)*v[0]+u*v[-1]
        return y

    def boosted_straight_initial_state(self,beta_cm_x=.25):
        if abs(beta_cm_x)>=1:
            raise ValueError('invalid boost')
        y=self.initial_state()
        x,v=self._unpack(y)
        b=self.params.initial_beta
        v[0,0]=(-b+beta_cm_x)/(1-b*beta_cm_x)
        v[-1,0]=(b+beta_cm_x)/(1+b*beta_cm_x)
        for i in range(1,self.n-1):
            u=i/(self.n-1)
            v[i]=(1-u)*v[0]+u*v[-1]
        x[:,0]*=math.sqrt(1-beta_cm_x**2)
        return y

    def constraints(self,state):
        x,v=self._unpack(state)
        dx=np.diff(x,axis=0)
        dv=np.diff(v,axis=0)
        ell=np.linalg.norm(dx,axis=1)
        if np.min(ell)<=20*self.params.regularizer_cm:
            raise MeshQualityViolation('equidistributed NG gauge reached near-zero segment')
        tangent=dx/ell[:,None]
        dellt=(C_CM_PER_PS*np.einsum('ij,ij->i',dv,tangent))
        g=ell[1:]-ell[:-1]
        gdot=dellt[1:]-dellt[:-1]
        # gddot/c = J*vdot + c*(|dv_perp|^2/ell differences).
        perp=dv-np.sum(dv*tangent,axis=1)[:,None]*tangent
        geom=C_CM_PER_PS*np.sum(perp*perp,axis=1)/ell
        b=-(geom[1:]-geom[:-1])
        J=np.zeros((self.n-2,3*self.n))
        for j in range(1,self.n-1):
            row=j-1
            J[row,3*(j-1):3*j]+=tangent[j-1]
            J[row,3*j:3*(j+1)]-=tangent[j-1]+tangent[j]
            J[row,3*(j+1):3*(j+2)]+=tangent[j]
        return J,b,g,gdot

    def rhs(self,t,state):
        n=self.n
        x,v=self._unpack(state)
        action=self.action_quantities(state,derivatives=True)
        M=action['mass_matrix']
        source=(action['forces_gev_per_ps']-
                action['convective_gev_per_ps']).ravel()
        J,b,g,gdot=self.constraints(state)
        k=n-2
        A=np.zeros((3*n+k,3*n+k))
        A[:3*n,:3*n]=M
        A[:3*n,3*n:]=J.T
        A[3*n:,:3*n]=J
        rhs=np.r_[source,b]
        try:
            solution=np.linalg.solve(A,rhs)
        except np.linalg.LinAlgError as exc:
            raise MeshQualityViolation(f'singular gauge-constrained NG KKT matrix: {exc}')
        vdot=solution[:3*n].reshape(n,3)
        if not np.all(np.isfinite(vdot)):
            raise MeshQualityViolation('nonfinite constrained NG acceleration')
        output=np.zeros(len(state))
        output[:3*n]=(C_CM_PER_PS*v).ravel()
        output[3*n:6*n]=vdot.ravel()
        output[-2]=math.sqrt(max(0,1-float(v[0]@v[0])))
        output[-1]=math.sqrt(max(0,1-float(v[-1]@v[-1])))
        return output

    def solve(self,state=None,duration_ps=None,strict_mesh=True,
              pinch_fraction_of_initial_length=4e-4,
              gauge_length_tol_cm=2e-5,
              gauge_length_rate_tol_cm_per_ps=2e-5):
        if state is None:
            state=self.initial_state()
        if duration_ps is None:
            duration_ps=self.params.duration_ps
        if strict_mesh:
            require_good_mesh(state,self.n,self.params.regularizer_cm)
        self.constraints(state)
        initial=self.action_quantities(state)
        # Mesh-independent PHYSICAL total arclength cutoff for uniform
        # segments. An element-specific cutoff with a fixed length would
        # falsely stop more-refined meshes sooner at the same encounter.
        # This is an engineering stop, NOT a confinement/annihilation scale.
        if not 0<pinch_fraction_of_initial_length<1:
            raise ValueError('pinch fraction must be in (0,1)')
        total_pinch_cm=(self.params.initial_length_cm*
                        pinch_fraction_of_initial_length)
        min_factor=max(100.,total_pinch_cm/
                       ((self.n-1)*self.params.regularizer_cm))
        events=(make_mesh_stop_events(self.n,self.params.regularizer_cm,
                                     minimum_factor=min_factor)
                if strict_mesh else None)
        sol=solve_ivp(self.rhs,(0,duration_ps),state,method='DOP853',
                      atol=1e-10,rtol=2e-9,max_step=self.params.max_step_ps,
                      events=events)
        if not sol.success:
            raise RuntimeError(sol.message)
        if strict_mesh:
            assert_no_mesh_stop(sol)
        final=self.action_quantities(sol.y[:,-1])
        J,b,g,gdot=self.constraints(sol.y[:,-1])
        max_g=float(np.max(np.abs(g)))
        max_gdot=float(np.max(np.abs(gdot)))
        if max_g>gauge_length_tol_cm or max_gdot>gauge_length_rate_tol_cm_per_ps:
            raise MeshQualityViolation(
                f'equal-arclength gauge drift: {max_g} cm and {max_gdot} cm/ps')
        return {
            'status':'EXPERIMENTAL_ZERO_INTERIOR_MASS_NG_GAUGE_ONLY',
            'duration_ps':duration_ps,'nodes':self.n,
            'initial_energy_gev':initial['energy_gev'],
            'final_energy_gev':final['energy_gev'],
            'energy_drift_gev':final['energy_gev']-initial['energy_gev'],
            'linear_momentum_drift_gev':
                float(np.linalg.norm(final['momentum_gev']-initial['momentum_gev'])),
            'angular_momentum_drift_hbar':
                float(np.linalg.norm(final['angular_momentum_over_hbar']-
                                     initial['angular_momentum_over_hbar'])),
            'max_equal_arclength_deviation_cm':max_g,
            'max_equal_arclength_rate_residual_cm_per_ps':max_gdot,
            'mesh_quality':self.mesh_quality(sol.y),
            'n_rhs_evaluations':sol.nfev,
            'limitations':[
                'midpoint finite element semidiscrete gauge, not continuum NG validated',
                'no physical remeshing or reconnection at all-endpoint crossing',
                'no weak recoil transport validated under new gauge',
                'same-time straight-string boost is an approximate kinematic slice'
            ]
        },sol
