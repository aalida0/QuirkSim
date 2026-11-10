"""Variationally discretized Nambu--Goto string with massive endpoints.

The target-time-gauge ACTION is discretized into straight finite elements:
  L = -sum_end m sqrt(1-v_end^2)
      -sum_interior mu_g sqrt(1-v_i^2)
      -T sum_segments |dx| sqrt(1-v_perp,mid^2),
  T = Lambda_IC^2 / (hbar*c) [GeV/cm].

An explicit tiny interior mass mu_g regularizes the tangent-direction
reparameterization zero modes. It is NOT physical and must be extrapolated
towards zero with spatial refinement. Segment midpoint quadrature and
piecewise-linear embedding do NOT yield an exact continuum NG solution.
Finite-element motion carries its own canonical energy/momentum/L.

Important: currently this module evolves undecayed endpoints only. The
weak-decay impulse coupling is deliberately BLOCKED pending full worldsheet
energy/momentum continuity validation. No claim about annihilation.
"""
from dataclasses import dataclass, replace
import math
import numpy as np
from scipy.integrate import solve_ivp
from .physics import benchmark
from .decays import sample_decay
from .dynamics import C_CM_PER_PS, HBARC_GEV_CM
from .field_dynamics import MAGNETIC_FORCE_GEV_PER_PS_T

IDENTITY=np.eye(3)


@dataclass(frozen=True)
class StringParameters:
    mass_gev: float = 225.0
    lambda_ev: float = 774.0
    nodes: int = 5
    initial_length_cm: float = 0.1
    initial_beta: float = .4
    gauge_mass_fraction: float = 1e-3
    regularizer_cm: float = 1e-9
    b_tesla: float = 0.0
    duration_ps: float = 3.0
    max_step_ps: float = 0.1

    def __post_init__(self):
        if self.mass_gev<=0 or self.lambda_ev<=0 or self.initial_length_cm<=0:
            raise ValueError('mass, confinement scale and length must be positive')
        if self.nodes<3 or self.nodes>32:
            raise ValueError('need 3..32 string vertices')
        if not 0<=self.initial_beta<1:
            raise ValueError('invalid endpoint beta')
        if self.gauge_mass_fraction<=0 or self.regularizer_cm<=0:
            raise ValueError('finite gauge regulator required')
        if self.duration_ps<=0 or self.max_step_ps<=0 or not math.isfinite(self.b_tesla):
            raise ValueError('invalid integration or field parameters')


class VariationalString:
    """Semidiscrete Euler-Lagrange solver: d/dt(dL/dv)=c*dL/dx+q v x B."""

    def __init__(self,params):
        self.params=params
        self.n=params.nodes
        self.tension=(params.lambda_ev*1e-9)**2/HBARC_GEV_CM
        m=benchmark('triplet',params.mass_gev)['charged_mass_gev']
        # This gauge mass must be extrapolated; never confuse with string density.
        self.gauge_mass=(params.gauge_mass_fraction*self.tension*
                         params.initial_length_cm/(self.n-2))
        self.masses=np.full(self.n,self.gauge_mass)
        self.masses[0]=self.masses[-1]=m
        self.charge=np.zeros(self.n)
        self.charge[0]=+1
        self.charge[-1]=-1

    def initial_state(self,transverse_velocity=0.0, middle_offset_cm=0.0):
        """A stretched-string slice, NOT the pointlike production time."""
        if abs(transverse_velocity)>=math.sqrt(1-self.params.initial_beta**2):
            raise ValueError('initial velocity exceeds c')
        n=self.n
        x=np.zeros((n,3))
        x[:,0]=np.linspace(-self.params.initial_length_cm/2,
                            self.params.initial_length_cm/2,n)
        if middle_offset_cm:
            x[n//2,1]=middle_offset_cm
        v=np.zeros((n,3))
        v[0,0]=-self.params.initial_beta
        v[-1,0]=self.params.initial_beta
        v[0,1]=transverse_velocity
        return np.r_[x.ravel(),v.ravel(),[0.,0.]]

    def _unpack(self,state):
        size=self.n*3
        return (np.asarray(state[:size]).reshape(self.n,3),
                np.asarray(state[size:2*size]).reshape(self.n,3))

    def action_quantities(self,state,derivatives=False):
        """Return conserved quantities and optional semidiscrete derivatives.

        Variation of the midpoint Nambu-Goto segment action gives canonical
        momenta including string degrees and the consistent mass matrix.
        """
        x,v=self._unpack(state)
        n=self.n
        speed2=np.einsum('ij,ij->i',v,v)
        if np.max(speed2)>=1:
            raise ValueError('superluminal discretized string vertex')
        gamma=1/np.sqrt(1-speed2)
        momenta=self.masses[:,None]*gamma[:,None]*v
        energy=float(np.sum(self.masses*gamma))
        M=np.zeros((3*n,3*n)) if derivatives else None
        forces=np.zeros((n,3)) if derivatives else None
        convective=np.zeros((n,3)) if derivatives else None
        if derivatives:
            for i in range(n):
                block=self.masses[i]*gamma[i]*(IDENTITY+gamma[i]**2*np.outer(v[i],v[i]))
                sl=slice(3*i,3*i+3)
                M[sl,sl]+=block
        for i in range(n-1):
            d=x[i+1]-x[i]
            ell=math.sqrt(float(d@d)+self.params.regularizer_cm**2)
            tangent=d/ell
            velocity=.5*(v[i]+v[i+1])
            upar=float(velocity@tangent)
            perp=velocity-upar*tangent
            h2=1-float(velocity@velocity)+upar**2
            if h2<=0:
                raise ValueError('invalid worldsheet segment time-likeness')
            h=math.sqrt(h2)
            T=self.tension
            shared=T*ell/(2*h)*perp
            momenta[i]+=shared
            momenta[i+1]+=shared
            energy+=T*ell/h
            if not derivatives:
                continue
            left=slice(3*i,3*i+3)
            right=slice(3*(i+1),3*(i+1)+3)
            Q=IDENTITY-np.outer(tangent,tangent)
            block=T*ell/4*(Q/h+np.outer(perp,perp)/h**3)
            for s in (left,right):
                for r in (left,right):
                    M[s,r]+=block
            force=C_CM_PER_PS*T*(h*tangent+(upar/h)*perp)
            forces[i]+=force
            forces[i+1]-=force
            # Change of canonical momentum due to moving x at FIXED v.
            d_dot=C_CM_PER_PS*(v[i+1]-v[i])
            ell_dot=float(tangent@d_dot)
            tan_dot=(d_dot-tangent*ell_dot)/ell
            upar_dot=float(velocity@tan_dot)
            perp_dot=-upar_dot*tangent-upar*tan_dot
            h_dot=upar*upar_dot/h
            partial=T/2*((ell_dot/h-ell*h_dot/h**2)*perp+
                         ell/h*perp_dot)
            convective[i]+=partial
            convective[i+1]+=partial
        total_P=np.sum(momenta,axis=0)
        total_J=np.sum(np.cross(x,momenta),axis=0)/HBARC_GEV_CM
        result={'energy_gev':energy,'momentum_gev':total_P,
                'angular_momentum_over_hbar':total_J,
                'canonical_vertex_momenta_gev':momenta}
        if derivatives:
            B=np.array([0.,0.,self.params.b_tesla])
            for i in (0,n-1):
                forces[i]+=self.charge[i]*MAGNETIC_FORCE_GEV_PER_PS_T*np.cross(v[i],B)
            result.update({'mass_matrix':M,'forces_gev_per_ps':forces,
                           'convective_gev_per_ps':convective})
        return result

    def rhs(self,t,state):
        n=self.n
        x,v=self._unpack(state)
        q=self.action_quantities(state,derivatives=True)
        vdot=np.linalg.solve(q['mass_matrix'],
                             (q['forces_gev_per_ps']-
                              q['convective_gev_per_ps']).ravel()).reshape(n,3)
        if not np.all(np.isfinite(vdot)):
            raise FloatingPointError('nonfinite acceleration')
        out=np.zeros(len(state))
        out[:3*n]=(C_CM_PER_PS*v).ravel()
        out[3*n:6*n]=vdot.ravel()
        out[-2]=math.sqrt(1-float(v[0]@v[0]))
        out[-1]=math.sqrt(1-float(v[-1]@v[-1]))
        return out

    def mesh_quality(self,states):
        """Heuristic spatial mesh-health *gate*, not a continuum guarantee.

        Reject under-resolved elements near the short-distance regulator
        and extreme adjacent material-vertex spacing. A conservative
        remeshing/worldsheet gauge algorithm is not yet implemented.
        """
        arr=np.asarray(states,float)
        dim=6*self.n+2
        if arr.ndim==1:
            if arr.shape[0]!=dim:
                raise ValueError('invalid state dimension')
            arr=arr[None,:]
        elif arr.ndim==2:
            if arr.shape[0]==dim:
                arr=arr.T
            elif arr.shape[1]!=dim:
                raise ValueError('invalid state dimension')
        else:
            raise ValueError('expected state or matrix of states')
        x=arr[:,:3*self.n].reshape(-1,self.n,3)
        length=np.linalg.norm(np.diff(x,axis=1),axis=2)
        min_l=float(length.min())
        worst=float(np.max(length.max(axis=1)/
                           np.maximum(length.min(axis=1),self.params.regularizer_cm)))
        short=min_l<100*self.params.regularizer_cm
        unbalanced=worst>50
        return {'minimum_segment_cm':min_l,
                'worst_segment_length_ratio':worst,
                'short_relative_to_regulator':bool(short),
                'length_ratio_above_50':bool(unbalanced),
                'quality_pass':bool(not short and not unbalanced),
                'status':'HEURISTIC_MESH_QUALITY_ONLY'}

    def weak_decay_impulse(self,state,endpoint,channel,seed=1,energy_tolerance_gev=1e-4,
                           strict=True):
        """Apply a weak endpoint transition by a CANONICAL momentum impulse.

        P_owning_vertex(after) = P_vertex(before) - sum(p_emitted).
        Other string vertices retain their canonical momenta; invert the
        discretized velocity->momentum relation to obtain updated velocities.

        The emitted products have exact on-shell four-momentum using the
        local parent *particle* momentum. The finite-element string can
        experience an O(discretization) discontinuity in total energy:
        the residual is explicitly calculated, never hidden. With
        strict=True, any residual > tolerance REJECTS and reverts the
        event. This is a development gate, not certified physical coupling.
        """
        if endpoint not in (0,self.n-1):
            raise ValueError('weak decays allowed only at string endpoints')
        if channel not in ('pion','electron_neutrino','muon_neutrino'):
            raise ValueError('unknown channel')
        if energy_tolerance_gev<0:
            raise ValueError('invalid energy tolerance')
        y=np.array(state,dtype=float,copy=True)
        x,v=self._unpack(y)
        parent_mass=float(self.masses[endpoint])
        daughter_mass=float(self.params.mass_gev)
        if parent_mass<=daughter_mass:
            raise ValueError('endpoint already neutral or inconsistent masses')
        initial=self.action_quantities(y)
        parent_particle_p=parent_mass*v[endpoint]/math.sqrt(1-float(v[endpoint]@v[endpoint]))
        outgoing=sample_decay(channel,parent_particle_p,parent_mass,daughter_mass,
                              np.random.default_rng(seed))
        emitted=[vec for name,vec in outgoing.items() if name!='daughter']
        p_emit=sum((vec[1:] for vec in emitted),np.zeros(3))
        e_emit=sum(float(vec[0]) for vec in emitted)
        p_target=initial['canonical_vertex_momenta_gev'].copy()
        p_target[endpoint]-=p_emit
        old_charge=float(self.charge[endpoint])
        self.masses[endpoint]=daughter_mass
        self.charge[endpoint]=0.
        v_guess=v.copy()
        daughter=outgoing['daughter']
        v_guess[endpoint]=daughter[1:]/daughter[0]
        sz=3*self.n
        try:
            for step in range(24):
                y[sz:2*sz]=v_guess.ravel()
                current=self.action_quantities(y,derivatives=True)
                diff=(p_target-current['canonical_vertex_momenta_gev']).ravel()
                if float(np.max(np.abs(diff)))<3e-10:
                    break
                change=np.linalg.solve(current['mass_matrix'],diff).reshape(self.n,3)
                factor=1.
                while factor>1e-8 and np.max(np.einsum('ij,ij->i',
                                    v_guess+factor*change,v_guess+factor*change))>=.999999:
                    factor*=.5
                if factor<=1e-8:
                    raise RuntimeError('momentum inversion leaves physical velocity domain')
                v_guess+=factor*change
            else:
                raise RuntimeError('finite-element momentum inversion failed to converge')
            y[sz:2*sz]=v_guess.ravel()
            after=self.action_quantities(y)
            mismatch=float(initial['energy_gev']-after['energy_gev']-e_emit)
            delta_P=np.array(initial['momentum_gev']-
                              after['momentum_gev']-p_emit)
            delta_J=np.array(initial['angular_momentum_over_hbar']-
                        after['angular_momentum_over_hbar']-
                        np.cross(x[endpoint],p_emit)/HBARC_GEV_CM)
            passed=abs(mismatch)<=energy_tolerance_gev
            if strict and not passed:
                raise RuntimeError('worldsheet weak-decay energy mismatch '
                                   f'{mismatch:.6g} GeV exceeds '
                                   f'{energy_tolerance_gev:.6g} GeV')
            diagnostics={
                'channel':channel,'endpoint':endpoint,'energy_residual_gev':mismatch,
                'emitted_energy_gev':e_emit,
                'momentum_residual_gev':delta_P.tolist(),
                'angular_momentum_residual_hbar':delta_J.tolist(),
                'energy_gate_pass':passed,
                'caveat':'finite-element instantaneous endpoint impulse; '
                          'worldsheet-limit and causal propagation not validated'
            }
            return y,diagnostics
        except Exception:
            self.masses[endpoint]=parent_mass
            self.charge[endpoint]=old_charge
            raise

    def solve(self,state=None, duration_ps=None,store_trace=False):
        if state is None:
            state=self.initial_state()
        if duration_ps is None:
            duration_ps=self.params.duration_ps
        initial=self.action_quantities(state)
        t_eval=np.linspace(0,duration_ps,41) if store_trace else None
        sol=solve_ivp(self.rhs,(0,duration_ps),state,method='DOP853',
                      atol=1e-10,rtol=2e-9,max_step=self.params.max_step_ps,
                      t_eval=t_eval)
        if not sol.success:
            raise RuntimeError(sol.message)
        final=self.action_quantities(sol.y[:,-1])
        diagnostics={
            'model':'variational_midpoint_Nambu_Goto_finite_elements',
            'vertices':self.n,'duration_ps':duration_ps,
            'gauge_mass_gev':self.gauge_mass,
            'initial_energy_gev':initial['energy_gev'],
            'final_energy_gev':final['energy_gev'],
            'energy_drift_gev':final['energy_gev']-initial['energy_gev'],
            'initial_momentum_gev':initial['momentum_gev'].tolist(),
            'final_momentum_gev':final['momentum_gev'].tolist(),
            'initial_angular_momentum_hbar':initial['angular_momentum_over_hbar'].tolist(),
            'final_angular_momentum_hbar':final['angular_momentum_over_hbar'].tolist(),
            'final_endpoint_positions_cm':self._unpack(sol.y[:,-1])[0][[0,-1]].tolist(),
            'n_rhs_evaluations':sol.nfev,
            'mesh_quality':self.mesh_quality(sol.y),
            'physics_gate':'SHORT_TIME_NUMERICAL_ONLY' if self.mesh_quality(sol.y)['quality_pass']
                           else 'BLOCKED_MESH_QUALITY',
            'limitations':[
                'midpoint finite-element Nambu-Goto action; not continuum validated',
                'interior tangent gauge regulator adds artificial massive degrees',
                'initialized at finite string length; not yet pointlike creation',
                'charged weak-decay impulse matching NOT implemented',
                'no radiation, annihilation, or detector materials']
        }
        return diagnostics,sol
