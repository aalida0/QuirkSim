"""Initialize triplet string evolution from actual Pythia8 endpoint four-vectors.

The original field_dynamics.simulate() assumes TWO CHARGED endpoints and a
fixed common boost. This independent entry point accepts the true charged-
neutral or charged-charged produced topology. Quirk charge is retained when
the string connects neutral endpoints. The model is still a straight-string
endpoint-force approximation, not a full worldsheet simulation.
"""
from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from .decays import sample_decay
from .dynamics import C_CM_PER_PS, sample_plan
from .field_dynamics import MAGNETIC_FORCE_GEV_PER_PS_T, string_force_pair
from .physics import benchmark

CHARGES = {9000111:1, -9000111:-1, 9000112:0, -9000112:0,
           9000113:-1, -9000113:1}
# Approximate ATLAS barrel pixel layer cylinder radii, in cm.
# NOT a sensor geometry: no finite z, endcaps, modules or hit efficiencies.
ATLAS_PIXEL_BARREL_RADII_CM = (3.325, 5.05, 8.85, 12.25)


@dataclass(frozen=True)
class ProductionParameters:
    lambda_ev: float = 774.
    b_tesla: float = 2.
    duration_ps: float = 200.
    max_step_ps: float = 2.
    regulator_cm: float = 1e-7
    seed: int = 10
    mode: str = 'recoil'
    force_model: str = 'ng_endpoint'
    mass_tolerance_gev: float = 0.002
    layer_radii_cm: tuple = ATLAS_PIXEL_BARREL_RADII_CM

    def __post_init__(self):
        if self.mode not in ('stable','no_recoil','recoil'):
            raise ValueError('invalid weak-decay mode')
        if self.force_model not in ('central_proxy','ng_endpoint'):
            raise ValueError('invalid string-force model')
        if min(self.lambda_ev,self.duration_ps,self.max_step_ps,
               self.regulator_cm,self.mass_tolerance_gev)<=0:
            raise ValueError('positive numerical and physical parameters required')
        if not math.isfinite(self.b_tesla):
            raise ValueError('invalid magnetic field')
        radii = tuple(self.layer_radii_cm)
        if (not radii or any(not math.isfinite(r) or r<=0 for r in radii)
                or any(a>=b for a,b in zip(radii,radii[1:]))):
            raise ValueError('layer radii must be positive, strictly increasing')



def read_events(path, limit=0):
    if limit<0:
        raise ValueError('limit must be >= 0')
    with Path(path).open(encoding='utf-8') as f:
        for index,line in enumerate(f,1):
            if not line.strip():
                continue
            try:
                event=json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f'Invalid JSONL event on line {index}') from exc
            yield event
            if limit and index>=limit:
                break


def validated_endpoints(event, expected_neutral_mass=225., tolerance=0.002):
    if event.get('schema_version')!=1:
        raise ValueError('Unsupported endpoint record schema')
    if event.get('stage')!='PYTHIA8_LO_QCD_ISR_FSR':
        raise ValueError('Expected Pythia8 showered endpoints')
    endpoints=event.get('endpoints',[])
    if len(endpoints)!=2:
        raise ValueError('Need precisely two produced quirk endpoints')
    if {int(x.get('quirk_number',0)) for x in endpoints}!={-1,1}:
        raise ValueError('The pair must consist of quirk and antiquark')
    bench=benchmark('triplet',expected_neutral_mass)
    for endpoint in endpoints:
        pdg=int(endpoint.get('pdg_id',0))
        if pdg not in CHARGES:
            raise ValueError('Unexpected endpoint PDG')
        if int(endpoint.get('charge_e',100))!=CHARGES[pdg]:
            raise ValueError('Electric charge mismatches PDG')
        if int(endpoint.get('quirk_number',0))!=(1 if pdg>0 else -1):
            raise ValueError('Quirk number mismatches PDG')
        p4=np.asarray(endpoint.get('p4_gev',[]),dtype=float)
        if p4.shape!=(4,) or not np.isfinite(p4).all():
            raise ValueError('Invalid endpoint four-momentum')
        mass=float(endpoint.get('mass_gev',0))
        expected=(bench['charged_mass_gev'] if CHARGES[pdg] else expected_neutral_mass)
        if abs(mass-expected)>tolerance:
            raise ValueError(f'UFO/Pythia mass {mass:.6g} differs from QuirkSim {expected:.6g}')
        if not math.isclose(p4[0],math.hypot(mass,float(np.linalg.norm(p4[1:]))),
                            rel_tol=0,abs_tol=0.01):
            raise ValueError('Endpoint four-vector is not on shell')
    return endpoints,bench


def simulate_produced(event, params: ProductionParameters):
    """Evolve a physically produced pair from one pointlike pp vertex.

    Electroweak pion/lepton weak decays occur ONLY for endpoints initially
    charged, with proper-time clocks. Quirk anti-particle sign is retained.
    The radiative recoil and hadron/matter response are not simulated.
    """
    endpoints, bench=validated_endpoints(event,225.,params.mass_tolerance_gev)
    charges=[int(x['charge_e']) for x in endpoints]
    charged=[q!=0 for q in charges]
    initially_charged=sum(charged)
    masses=[float(x['mass_gev']) for x in endpoints]
    y=np.zeros(14,dtype=float)
    for i,end in enumerate(endpoints):
        y[6+3*i:9+3*i]=np.asarray(end['p4_gev'][1:],float)
    plan=sample_plan(bench['lifetime_ps'],bench['widths_gev'],params.seed)
    tension=(params.lambda_ev*1e-9)**2
    t=0.
    decays=[]
    minima=[]
    # The first outward crossing is recorded even if it occurs *after* decay.
    # A charged crossing is not a detector hit: barrel length is unbounded here.
    first_outward_crossings = {}
    first_decay_time=None
    bvec=np.array([0.,0.,params.b_tesla])
    while t < params.duration_ps-1e-10:
        def rhs(_t,z):
            out=np.zeros(14)
            force1,force2=string_force_pair(z,masses,tension,
                 params.regulator_cm,params.force_model)
            for i,force in enumerate((force1,force2)):
                pi=z[6+3*i:9+3*i]
                energy=math.hypot(masses[i],float(np.linalg.norm(pi)))
                vel=pi/energy
                out[3*i:3*i+3]=C_CM_PER_PS*vel
                out[6+3*i:9+3*i]=force+charges[i]*MAGNETIC_FORCE_GEV_PER_PS_T*np.cross(vel,bvec)
                out[12+i]=masses[i]/energy
            return out

        def close_min(_t,z):
            r=z[:3]-z[3:6]
            v0=z[6:9]/math.hypot(masses[0],float(np.linalg.norm(z[6:9])))
            v1=z[9:12]/math.hypot(masses[1],float(np.linalg.norm(z[9:12])))
            return float(np.dot(r,v0-v1))
        close_min.direction=1
        close_min.terminal=False
        active=[i for i in range(2) if charges[i]!=0 and params.mode!='stable']
        def make_decay_clock(i):
            def clock(_t,z):
                return z[12+i]-plan[i]['proper_time_ps']
            clock.direction=1
            clock.terminal=True
            return clock

        # Integrator root finds: outward crossings of infinite barrel
        # cylinders, for both endpoints whether or not they are charged.
        layer_keys=[(i,float(radius)) for i in range(2)
                    for radius in params.layer_radii_cm]
        def make_radius_crossing(i,radius):
            def cross(_t,z):
                return math.hypot(z[3*i],z[3*i+1])-radius
            cross.direction=1
            cross.terminal=False
            return cross

        sol=solve_ivp(rhs,(t,params.duration_ps),y,method='DOP853',
            events=([close_min]+[make_decay_clock(i) for i in active]
                    +[make_radius_crossing(i,r) for i,r in layer_keys]),
            rtol=2e-9,atol=1e-10,max_step=params.max_step_ps)
        if not sol.success:
            raise RuntimeError(sol.message)
        for tt,yy in zip(sol.t_events[0],sol.y_events[0]):
            if tt>1e-6:
                minima.append({'time_ps':float(tt),'separation_cm':
                              float(np.linalg.norm(yy[:3]-yy[3:6]))})
        first_layer_idx=1+len(active)
        for j,(i,radius) in enumerate(layer_keys):
            key=(i,radius)
            if key in first_outward_crossings:
                continue
            if not len(sol.t_events[first_layer_idx+j]):
                continue
            tt=float(sol.t_events[first_layer_idx+j][0])
            zz=sol.y_events[first_layer_idx+j][0]
            first_outward_crossings[key]={
                'endpoint':i, 'radius_cm':radius, 'time_ps':tt,
                'z_cm':float(zz[3*i+2]),
                'phi_rad':math.atan2(zz[3*i+1],zz[3*i]),
                'charge_at_crossing_e':int(charges[i]),
                'initial_charge_e':int(endpoints[i]['charge_e']),
            }
        y=sol.y[:,-1].copy()
        t=float(sol.t[-1])
        triggered=[idx for idx in range(len(active)) if len(sol.t_events[idx+1])]
        if not triggered:
            break
        i=active[triggered[0]]
        channel=plan[i]['channel']
        # The entire charged->neutral recoil is applied to the endpoint;
        # emitted weak products are recorded but not tracked through material.
        products=sample_decay(channel,y[6+3*i:9+3*i],masses[i],225.,
                              np.random.default_rng(plan[i]['rest_seed']))
        if params.mode=='recoil':
            y[6+3*i:9+3*i]=products['daughter'][1:]
        emitted={key:list(map(float,vec)) for key,vec in products.items()
                 if key!='daughter'}
        decays.append({'endpoint':i,'parent_pdg':int(endpoints[i]['pdg_id']),
                       'emitted_charge_e':charges[i],
                       'time_ps':t,'channel':channel,
                       'position_cm':[float(c) for c in y[3*i:3*i+3]],
                       'radius_cm':math.hypot(y[3*i],y[3*i+1]),
                       'emitted_p4_gev':emitted})
        charges[i]=0
        masses[i]=225.
        if first_decay_time is None: first_decay_time=t

    later=[x for x in minima if first_decay_time is not None
           and x['time_ps']>first_decay_time+1e-7]
    return {
        'schema_version':1,'production_event_index':event.get('event_index'),
        'initial_pdg_ids':[int(e['pdg_id']) for e in endpoints],
        'initial_p4_gev':[list(map(float,e['p4_gev'])) for e in endpoints],
        'initial_charges_e':[int(e['charge_e']) for e in endpoints],
        'initial_charged_endpoints':initially_charged,
        'final_charged_endpoints':sum(q!=0 for q in charges),
        'n_weak_decays':len(decays),'decays':decays,'turning_points':minima,
        'barrel_layer_radii_cm':list(params.layer_radii_cm),
        'first_outward_barrel_crossings':
            sorted(first_outward_crossings.values(),
                   key=lambda x:(x['radius_cm'],x['endpoint'])),
        'post_decay_closest_separation_cm':
            min((x['separation_cm'] for x in later),default=None),
        'end_time_ps':t,'pair_pt_gev':event['pair_pt_gev'],
        'final_endpoint_p3_gev':[list(map(float,y[6:9])),list(map(float,y[9:12]))],
        'params':vars(params),
        'caveats':['Pythia8 QCD ISR/FSR only, no infracolor emission',
                   'straight-string endpoint-force approximation, not worldsheet',
                   'no annihilation, energy loss or detector response',
                   'initially neutral endpoints cannot decay again',
                   'infinite barrel cylinders, no z acceptance, endcaps or sensor hits',
                   'no per-layer charge deposition or reconstruction efficiency']
    }
