"""Event-driven two-endpoint *central string proxy* with complete LO weak decays.

UNITS: GeV mass/momentum/energy, cm positions, ps lab time. The potential
V(r) = sigma [sqrt(r^2+a^2)-a] has correct long-string tension, but is
NOT a moving Nambu--Goto string worldsheet. Radiation/annihilation are absent.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.integrate import solve_ivp
from .physics import benchmark, HBAR_GEV_PS
from .decays import sample_decay

C_CM_PER_PS = 0.0299792458
HBARC_GEV_CM = HBAR_GEV_PS * C_CM_PER_PS


@dataclass(frozen=True)
class Parameters:
    representation: str = 'triplet'
    mass_neutral_gev: float = 225.0
    lambda_ev: float = 774.0
    beta: float = 0.4
    impact_cm: float = 0.0
    regulator_cm: float = 1e-7
    duration_ps: float = 500.0
    seed: int = 1
    mode: str = 'recoil'
    store_trace: bool = False
    max_step_ps: float = 2.0

    def __post_init__(self):
        if self.representation not in ('triplet', 'doublet'):
            raise ValueError('unknown electroweak representation')
        if self.mode not in ('stable', 'no_recoil', 'recoil'):
            raise ValueError('unknown simulation mode')
        if self.mass_neutral_gev <= 0 or self.lambda_ev <= 0 or self.duration_ps <= 0:
            raise ValueError('positive mass, Lambda and duration required')
        if not 0 < self.beta < 1 or self.impact_cm < 0 or self.regulator_cm <= 0 or self.max_step_ps <= 0:
            raise ValueError('invalid initial conditions / integration settings')


def sample_plan(lifetime, widths, seed):
    """Shared latent decay clock, channel and rest-frame random seed per endpoint.

    Passing the same seed in no_recoil/recoil guarantees same proper-time
    thresholds, channel assignment and rest-frame product kinematics.
    Lab decay times can differ physically because the endpoint velocities differ.
    """
    rng = np.random.default_rng(seed)
    channels = list(widths)
    fractions = np.array([widths[c] for c in channels], float)
    fractions = fractions / fractions.sum()
    tau = rng.exponential(lifetime, 2)
    channel_ids = rng.choice(len(channels), 2, p=fractions)
    rest_seeds = rng.integers(0, 2**32, size=2)
    return [{'proper_time_ps':float(tau[i]),
             'channel':channels[int(channel_ids[i])],
             'rest_seed':int(rest_seeds[i])} for i in range(2)]


def total_energy(y, masses, sigma, regulator_cm):
    r = np.linalg.norm(y[:3]-y[3:6])
    vstring = sigma / HBARC_GEV_CM * (math.hypot(r,regulator_cm)-regulator_cm)
    return (math.hypot(masses[0],np.linalg.norm(y[6:9])) +
            math.hypot(masses[1],np.linalg.norm(y[9:12])) + vstring)


def pion_two_body(momentum, parent_mass, daughter_mass, rng):
    """Backward-compatible pion sampler; returns daughter momentum and pion 4-vector."""
    products = sample_decay('pion',momentum,parent_mass,daughter_mass,rng)
    return products['daughter'][1:], (float(products['pion'][0]),products['pion'][1:].copy())


def propagate(parameters):
    """Propagate both endpoints with full pion/e/muon transitions and encounter roots.

    'closest_post_decay_turning_point_cm' is the minimum over local minima of
    |x1-x2| *after* first weak decay. It is never the trivial production point.
    Contact criterion 1/Lambda is geometric; NOT an annihilation probability.
    """
    p = parameters
    b = benchmark(p.representation,p.mass_neutral_gev)
    neutral_m = p.mass_neutral_gev
    charged_m = b['charged_mass_gev']
    sigma = (p.lambda_ev*1e-9)**2
    plan = sample_plan(b['lifetime_ps'], b['widths_gev'],p.seed)
    masses = [charged_m,charged_m]
    charged = [True,True]
    y = np.zeros(14)  # pos1 xyz, pos2 xyz, mom1 xyz, mom2 xyz, tau1, tau2
    y[:3] = [p.impact_cm/2,0,5e-8]
    y[3:6] = [-p.impact_cm/2,0,-5e-8]
    p_abs = charged_m*p.beta/math.sqrt(1-p.beta**2)
    y[6:9] = [0,0,p_abs]
    y[9:12] = [0,0,-p_abs]
    t=0.
    event_log=[]
    turning_points=[]
    trace=[]
    energies=[]
    first_decay_time=None
    min_global=float(np.linalg.norm(y[:3]-y[3:6]))
    max_global=min_global

    while t < p.duration_ps-1e-10:
        initial_time=t
        def rhs(_t,z):
            out=np.zeros(14)
            r=z[:3]-z[3:6]
            rnorm=math.hypot(float(np.linalg.norm(r)),p.regulator_cm)
            dp1=-sigma/HBAR_GEV_PS*r/rnorm
            out[6:9]=dp1
            out[9:12]=-dp1
            for i in range(2):
                momentum=z[6+3*i:9+3*i]
                energy=math.hypot(masses[i],np.linalg.norm(momentum))
                out[3*i:3*i+3]=C_CM_PER_PS*momentum/energy
                out[12+i]=masses[i]/energy
            return out

        def minimum_event(_t,z):
            r=z[:3]-z[3:6]
            v1=C_CM_PER_PS*z[6:9]/math.hypot(masses[0],np.linalg.norm(z[6:9]))
            v2=C_CM_PER_PS*z[9:12]/math.hypot(masses[1],np.linalg.norm(z[9:12]))
            return float(r@(v1-v2))
        minimum_event.direction=1
        minimum_event.terminal=False

        active=[i for i in range(2) if charged[i] and p.mode!='stable']
        def decay_event_factory(i):
            def dec_ev(_t,z):
                return z[12+i]-plan[i]['proper_time_ps']
            dec_ev.direction=1
            dec_ev.terminal=True
            return dec_ev

        energy_start=total_energy(y,masses,sigma,p.regulator_cm)
        sol=solve_ivp(rhs,(t,p.duration_ps),y,method='DOP853',
                      events=[minimum_event]+[decay_event_factory(i) for i in active],
                      atol=1e-10,rtol=2e-9,max_step=p.max_step_ps)
        if not sol.success:
            raise RuntimeError(sol.message)
        t=float(sol.t[-1]); y=sol.y[:,-1].copy()
        e_end=total_energy(y,masses,sigma,p.regulator_cm)
        energies.extend([energy_start,e_end])
        if p.store_trace:
            trace.extend([[float(ti),float(np.linalg.norm(zi[:3]-zi[3:6]))]
                          for ti,zi in zip(sol.t,sol.y.T) if not trace or ti>trace[-1][0]])
        separations=np.linalg.norm(sol.y[:3,:]-sol.y[3:6,:],axis=0)
        min_global=min(min_global,float(separations.min()))
        max_global=max(max_global,float(separations.max()))
        for tt,zz in zip(sol.t_events[0],sol.y_events[0]):
            if tt < 1e-6:
                continue
            turning_points.append({'time_ps':float(tt),
                                   'separation_cm':float(np.linalg.norm(zz[:3]-zz[3:6]))})
        fired=[j for j in range(len(active)) if len(sol.t_events[j+1])]
        if not fired:
            break
        i=active[fired[0]]
        channel=plan[i]['channel']
        products=sample_decay(channel,y[6+3*i:9+3*i],masses[i],neutral_m,
                              np.random.default_rng(plan[i]['rest_seed']))
        before_energy=total_energy(y,masses,sigma,p.regulator_cm)
        if p.mode=='recoil':
            y[6+3*i:9+3*i]=products['daughter'][1:]
        masses[i]=neutral_m
        charged[i]=False
        after_energy=total_energy(y,masses,sigma,p.regulator_cm)
        emitted={key: [float(x) for x in vec]
                 for key,vec in products.items() if key != 'daughter'}
        emitted_energy=sum(vec[0] for vec in emitted.values())
        energy_residual=(before_energy-after_energy-emitted_energy) if p.mode=='recoil' else None
        if p.mode=='recoil' and abs(energy_residual)>2e-8:
            raise RuntimeError(f'full-system energy discontinuity: {energy_residual}')
        event_log.append({'endpoint':i,'time_ps':float(t),
                          'proper_time_ps':float(y[12+i]),'channel':channel,
                          'emitted_fourvectors_gev':emitted,
                          'recoil_applied':p.mode=='recoil',
                          'energy_balance_residual_gev':energy_residual})
        if first_decay_time is None:
            first_decay_time=t
        if t <= initial_time+1e-13:
            raise RuntimeError('non-progressing event loop')

    if p.store_trace and trace:
        assert abs(trace[-1][0]-t)<1e-6
    encounters=[pt for pt in turning_points if first_decay_time is not None and
                pt['time_ps'] > first_decay_time+1e-7]
    core_radius_cm=HBARC_GEV_CM/(p.lambda_ev*1e-9)
    closest=min((pt['separation_cm'] for pt in encounters),default=None)
    out={
        'model':p.representation,'mode':p.mode,'mass_neutral_gev':neutral_m,
        'lambda_ev':p.lambda_ev,'lifetime_ps':b['lifetime_ps'],
        'delta_mev':b['delta_mev'],'branching_fractions':
            {k:v/sum(b['widths_gev'].values()) for k,v in b['widths_gev'].items()},
        'plan':plan, 'events':event_log,'end_time_ps':t,'charged':charged,
        'segment_energies_gev':energies,
        'minimum_separation_including_production_cm':min_global,
        'maximum_separation_cm':max_global,
        'first_decay_time_ps':first_decay_time,
        'turning_points':turning_points,
        'post_decay_encounter_count':len(encounters),
        'closest_post_decay_turning_point_cm':closest,
        'core_radius_cm':core_radius_cm,
        'post_decay_geometric_core_encounters':sum(pt['separation_cm']<=core_radius_cm for pt in encounters),
        'caveats':['central-string proxy, not Nambu-Goto',
                   'no radiation or annihilation',
                   'no-recoil counterfactual is unphysical']
    }
    if p.store_trace:
        out['trace_time_ps_separation_cm']=trace
    return out
