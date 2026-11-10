"""Quirk endpoint motion in a solenoid, with an approximate moving-string force.

Physics sources:
- Li et al., EPJC 80 (2020) 651, arXiv:2002.07503, Eqs. (2.1)-(2.3).
- ATLAS central solenoid nominal field B_z ~ 2 Tesla.

The 'ng_endpoint' option implements the published straight-string Nambu-Goto
*endpoint-force approximation*. It DOES NOT evolve the string worldsheet or
its intrinsic energy/momentum. It is NOT a full Nambu-Goto PDE solution.

Units: cm positions, ps time, GeV momenta, B in Tesla, c=0.0299792458 cm/ps.
"""
from dataclasses import dataclass, replace
import math
import numpy as np
from scipy.integrate import solve_ivp

from .physics import benchmark, HBAR_GEV_PS
from .decays import sample_decay
from .dynamics import C_CM_PER_PS, HBARC_GEV_CM, sample_plan

# Lorentz dp/dt: 0.299792458 [GeV/(T*m)] * c[m/ps] * (v/c x B).
C_M_PER_PS = C_CM_PER_PS / 100.0
MAGNETIC_FORCE_GEV_PER_PS_T = 0.299792458 * C_M_PER_PS


@dataclass(frozen=True)
class FieldParameters:
    model: str = 'triplet'
    mass_gev: float = 225.0
    lambda_ev: float = 774.0
    beta_rel: float = 0.4
    beta_cm_x: float = 0.25
    b_tesla: float = 2.0
    force_model: str = 'ng_endpoint'
    duration_ps: float = 1000.0
    seed: int = 10
    mode: str = 'recoil'
    regulator_cm: float = 1e-7
    max_step_ps: float = 2.0

    def __post_init__(self):
        if self.model not in ('triplet','doublet'):
            raise ValueError('invalid representation')
        if self.force_model not in ('central_proxy','ng_endpoint'):
            raise ValueError('invalid string force')
        if self.mode not in ('stable','no_recoil','recoil'):
            raise ValueError('invalid weak-decay mode')
        if not (0<self.beta_rel<1 and abs(self.beta_cm_x)<1):
            raise ValueError('invalid velocities')
        if self.mass_gev<=0 or self.lambda_ev<=0 or self.duration_ps<=0:
            raise ValueError('nonpositive benchmark')
        if self.regulator_cm<=0 or self.max_step_ps<=0 or not math.isfinite(self.b_tesla):
            raise ValueError('invalid numerical or magnetic inputs')


def ng_force(velocity, outward_unit, sigma_over_hbar):
    """Li et al. (2020), Eq. (2.2), in GeV/ps.

    outward_unit points from the opposite endpoint to this endpoint;
    force points inward, including transverse-velocity corrections.
    """
    vpar=float(velocity@outward_unit)
    vperp=velocity-vpar*outward_unit
    root=math.sqrt(max(1e-15,1-float(vperp@vperp)))
    return -sigma_over_hbar*(root*outward_unit+(vpar/root)*vperp)


def central_force(outward_unit, sigma_over_hbar):
    return -sigma_over_hbar*outward_unit


def initial_state(p, charged_mass):
    """Pointlike pair, exactly zero initial relative orbital angular momentum.

    In the pair CoM, momenta are +/-p* along x. Apply a common Lorentz
    boost beta_cm_x along x (a possible ISR-induced lab-frame pair boost).
    The charges are +1 and -1. Initial positions are a tiny regulated
    separation along the same x axis, NOT an arbitrary transverse offset.
    """
    y=np.zeros(14)
    y[0]=+p.regulator_cm/2
    y[3]=-p.regulator_cm/2
    pre=charged_mass*p.beta_rel/math.sqrt(1-p.beta_rel**2)
    en=charged_mass/math.sqrt(1-p.beta_rel**2)
    gam=1/math.sqrt(1-p.beta_cm_x**2)
    y[6]=gam*(pre+p.beta_cm_x*en)
    y[9]=gam*(-pre+p.beta_cm_x*en)
    return y


def relative_angular_momentum(y):
    """L_rel = (r1-r2) x (p1-p2)/2, in GeV*cm.

    A kinematic diagnostic. In a moving string the string itself can carry
    angular momentum: this is NOT full-system angular momentum.
    """
    r=y[:3]-y[3:6]
    relative_momentum=.5*(y[6:9]-y[9:12])
    return np.cross(r,relative_momentum)


def string_force_pair(y, masses, tension_gev2, regulator_cm, model):
    r=y[:3]-y[3:6]
    norm=float(np.linalg.norm(r))
    scale=norm/math.hypot(norm,regulator_cm)
    unit=r/norm if norm>0 else np.array([1.,0.,0.])
    coefficient=tension_gev2/HBAR_GEV_PS*scale
    if model=='central_proxy':
        return central_force(unit,coefficient),central_force(-unit,coefficient)
    v1=y[6:9]/math.hypot(masses[0],np.linalg.norm(y[6:9]))
    v2=y[9:12]/math.hypot(masses[1],np.linalg.norm(y[9:12]))
    return ng_force(v1,unit,coefficient),ng_force(v2,-unit,coefficient)


def simulate(p):
    """Run one pair with magnetic bending, both endpoint weak decays and encounters."""
    bench=benchmark(p.model,p.mass_gev)
    charged_mass=bench['charged_mass_gev']
    masses=[charged_mass,charged_mass]
    sigma=(p.lambda_ev*1e-9)**2
    charged=[True,True]
    plan=sample_plan(bench['lifetime_ps'],bench['widths_gev'],p.seed)
    y=initial_state(p,charged_mass)
    t=0.
    decays=[]
    minima=[]
    max_pre_L=0.
    max_post_L=0.
    first_decay_L=None
    energy_segments=[]

    def e_proxy(state):
        r=np.linalg.norm(state[:3]-state[3:6])
        return (math.hypot(masses[0],np.linalg.norm(state[6:9]))
                +math.hypot(masses[1],np.linalg.norm(state[9:12]))
                +sigma/HBARC_GEV_CM*(math.hypot(r,p.regulator_cm)-p.regulator_cm))

    while t < p.duration_ps-1e-10:
        pre_first=not decays

        def rhs(_t,state):
            out=np.zeros(14)
            f1,f2=string_force_pair(state,masses,sigma,p.regulator_cm,p.force_model)
            for i, f in enumerate((f1,f2)):
                pi=state[6+3*i:9+3*i]
                e=math.hypot(masses[i],np.linalg.norm(pi))
                v=pi/e
                out[3*i:3*i+3]=C_CM_PER_PS*v
                if charged[i]:
                    # B = +Bz; q1 = +1, q2 = -1, 0 after weak decay.
                    f=f+((1 if i==0 else -1)*MAGNETIC_FORCE_GEV_PER_PS_T*
                         np.cross(v,np.array([0.,0.,p.b_tesla])))
                out[6+3*i:9+3*i]=f
                out[12+i]=masses[i]/e
            return out

        def close_min(_t,state):
            r=state[:3]-state[3:6]
            v1=state[6:9]/math.hypot(masses[0],np.linalg.norm(state[6:9]))
            v2=state[9:12]/math.hypot(masses[1],np.linalg.norm(state[9:12]))
            return float(r@(v1-v2))
        close_min.direction=1
        close_min.terminal=False

        active=[i for i in range(2) if charged[i] and p.mode!='stable']
        def make_ev(i):
            def ev(_t,state):
                return state[12+i]-plan[i]['proper_time_ps']
            ev.direction=1
            ev.terminal=True
            return ev

        before=e_proxy(y)
        sol=solve_ivp(rhs,(t,p.duration_ps),y,method='DOP853',
            events=[close_min]+[make_ev(i) for i in active],
            atol=1e-10,rtol=2e-9,max_step=p.max_step_ps)
        if not sol.success:
            raise RuntimeError(sol.message)
        ang=np.cross(sol.y[:3,:].T-sol.y[3:6,:].T,
                     .5*(sol.y[6:9,:].T-sol.y[9:12,:].T))
        maximum=float(np.max(np.linalg.norm(ang,axis=1)))/HBARC_GEV_CM
        if pre_first:
            max_pre_L=max(max_pre_L,maximum)
        else:
            max_post_L=max(max_post_L,maximum)
        for tt,yy in zip(sol.t_events[0],sol.y_events[0]):
            if tt>1e-6:
                minima.append({'time_ps':float(tt),'separation_cm':
                               float(np.linalg.norm(yy[:3]-yy[3:6]))})
        y=sol.y[:,-1].copy()
        t=float(sol.t[-1])
        energy_segments.extend([before,e_proxy(y)])
        triggered=[k for k in range(len(active)) if len(sol.t_events[k+1])]
        if not triggered:
            break
        i=active[triggered[0]]
        oldL=relative_angular_momentum(y)/HBARC_GEV_CM
        if first_decay_L is None:
            first_decay_L=[float(z) for z in oldL]
        channel=plan[i]['channel']
        products=sample_decay(channel,y[6+3*i:9+3*i],masses[i],p.mass_gev,
                              np.random.default_rng(plan[i]['rest_seed']))
        if p.mode=='recoil':
            y[6+3*i:9+3*i]=products['daughter'][1:]
        charged[i]=False
        masses[i]=p.mass_gev
        updatedL=relative_angular_momentum(y)/HBARC_GEV_CM
        decays.append({'endpoint':i,'time_ps':t,'channel':channel,
                       'L_before_hbar':[float(z) for z in oldL],
                       'L_after_hbar':[float(z) for z in updatedL],
                       'delta_L_from_decay_hbar':
                           [float(z) for z in updatedL-oldL],
                       'recoil_applied':p.mode=='recoil'})

    first_time=decays[0]['time_ps'] if decays else None
    post=[r for r in minima if first_time is not None and r['time_ps']>first_time+1e-7]
    closest=min((r['separation_cm'] for r in post),default=None)
    core=HBARC_GEV_CM/(p.lambda_ev*1e-9)
    return {
      'model':p.model,'force_model':p.force_model,'mode':p.mode,
      'b_tesla':p.b_tesla,'beta_cm_x':p.beta_cm_x,'beta_rel':p.beta_rel,
      'seed':p.seed,'decays':decays,
      'plan':plan,
      'initial_relative_L_over_hbar':float(np.linalg.norm(relative_angular_momentum(initial_state(p,charged_mass))))/HBARC_GEV_CM,
      'max_pre_first_decay_L_over_hbar':max_pre_L,
      'first_decay_L_over_hbar':first_decay_L,
      'max_post_first_decay_L_over_hbar':max_post_L,
      'closest_post_decay_turning_point_cm':closest,
      'geometric_core_radius_cm':core,
      'post_decay_core_encounters':sum(r['separation_cm']<core for r in post),
      'turning_points':minima,
      'energy_proxy_segment_endpoints':energy_segments,
      'end_time_ps':t,
      'caveats':['straight-string Nambu-Goto endpoint FORCE approximation only',
                 'no worldsheet PDE or string energy-momentum degrees of freedom',
                 'no ionization, infracolor radiation or annihilation',
                 'uniform solenoid field; source at pointlike vertex',
                 'field response of emitted weak-decay products omitted']
    }
