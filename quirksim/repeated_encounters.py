"""Declared secular central-orbit model after a charge-changing handoff.

This is a phenomenological reduction, not NG string or quantum evolution.
E is pair excess energy (GeV); L is orbital angular momentum in hbar units.
Mean loss follows epsilon*Lambda per full signed radial reference period.
"""
from dataclasses import dataclass
import math
from scipy.optimize import brentq
from .physics import HBAR_GEV_PS
from .radiation_annihilation import OrbitAveragedModel,linear_cooling_time_ps


def pericenter_gev_inverse(mass,energy,angular,lambda_ev):
    """Inner turning point of H=2 sqrt(m^2+p_r^2+L^2/r^2)+sigma*r.

    Exact turning geometry in the specified central model. The secular clock
    below still uses the published radial period approximation.
    """
    if (not all(map(math.isfinite,[mass,energy,angular,lambda_ev]))
            or mass<=0 or energy<=0 or angular<0 or lambda_ev<=0):
        raise ValueError('positive mass/energy/scale and nonnegative L required')
    if angular==0:return 0.
    sigma=(lambda_ev*1e-9)**2;q=energy/(4*mass)
    a=angular*sigma/(energy*math.sqrt(mass*energy))
    def g(x):return x*x*(1-x)*(1+q*(1-x))-a*a
    peak=brentq(lambda x:4*q*x*x-3*(1+2*q)*x+2*(1+q),0.,1.)
    if g(peak)<0:raise ValueError('energy and angular momentum admit no central orbit')
    if g(peak)==0:return peak*energy/sigma
    return brentq(g,0.,peak,xtol=1e-300,rtol=1e-12)*energy/sigma


@dataclass(frozen=True)
class EncounterPrescription:
    mass_gev: float = 225.
    lambda_ev: float = 774.
    epsilon_per_full_period: float = .1
    gate_radius_multiple: float = 1.
    angular_loss_per_effective_emission: float = 0.
    capture_energy_in_lambda: float = 8.
    capture_max_L: float = 1.
    captured_binding_energy_gev: float = .050708354
    captured_singlet_fraction: float = .25
    singlet_annihilation_width_gev: float = 1.975004798e-5
    annihilation_optical_depth_per_eligible_encounter: float = 0.
    energy_step_fraction: float = .02

    def __post_init__(self):
        if not all(math.isfinite(x) for x in self.__dict__.values()):raise ValueError('finite model inputs required')
        if (self.mass_gev<=0 or self.lambda_ev<=0 or not 0<=self.epsilon_per_full_period<=1
                or self.gate_radius_multiple<=0 or not 0<=self.angular_loss_per_effective_emission<=self.gate_radius_multiple
                or self.capture_energy_in_lambda<=0 or self.capture_max_L<0
                or self.captured_binding_energy_gev<0 or not 0<=self.captured_singlet_fraction<=1
                or self.singlet_annihilation_width_gev<=0 or self.annihilation_optical_depth_per_eligible_encounter<0
                or not 0<self.energy_step_fraction<=.1):
            raise ValueError('invalid secular/capture prescription')


def evolve(initial_energy_gev,initial_L,horizon_ps,model):
    if not all(math.isfinite(x) for x in (initial_energy_gev,initial_L,horizon_ps)) or initial_energy_gev<=0 or initial_L<0 or horizon_ps<0:
        raise ValueError('physical initial state and nonnegative horizon required')
    lam=model.lambda_ev*1e-9;gate=model.gate_radius_multiple/lam
    cap=model.capture_energy_in_lambda*lam
    clock=OrbitAveragedModel(model.mass_gev,model.lambda_ev,model.epsilon_per_full_period,
        lam,0.,0.,1.,'Li2024 Eq7 mean loss; declared full signed radial clock',charge_magnitude=0.,include_photons=False)
    e=initial_energy_gev;l=initial_L;t=0.;radiated=0.;removed_L=0.;encounters=0.;n=0.
    history=[];status='HORIZON_REACHED';capture_time=None
    while True:
        r=pericenter_gev_inverse(model.mass_gev,e,l,model.lambda_ev)
        history.append({'time_ps':t,'energy_gev':e,'L':l,'pericenter_gev_inverse':r})
        if e<=cap*(1+1e-12):
            if l<=model.capture_max_L:
                status='CAPTURED';capture_time=t
            else:status='CAPTURE_ANGULAR_THRESHOLD_NOT_MET'
            break
        if t>=horizon_ps:break
        if r>gate:
            status='EMISSION_GATE_CLOSED';break
        if model.epsilon_per_full_period==0:
            status='NO_RADIATION';break
        target=max(cap,e*(1-model.energy_step_fraction))
        def next_L(energy):return max(0.,l-model.angular_loss_per_effective_emission*(e-energy)/lam)
        try:new_r=pericenter_gev_inverse(model.mass_gev,target,next_L(target),model.lambda_ev)
        except ValueError:
            status='SECULAR_UPDATE_OUTSIDE_ALLOWED_CENTRAL_ORBITS';break
        gate_closes=new_r>gate
        if gate_closes:
            target=brentq(lambda energy:pericenter_gev_inverse(model.mass_gev,energy,next_L(energy),model.lambda_ev)-gate,
                          target,e,xtol=1e-16,rtol=1e-12)
        dt=linear_cooling_time_ps(clock,e,target)
        if t+dt>horizon_ps:
            target=brentq(lambda energy:linear_cooling_time_ps(clock,e,energy)-(horizon_ps-t),
                          target,e*(1-1e-15),xtol=1e-16,rtol=1e-12)
            dt=horizon_ps-t;gate_closes=False
        new_l=next_L(target);count=(e-target)/lam
        # The reference full signed period has two radial encounters.
        encounters+=2*count/model.epsilon_per_full_period
        radiated+=e-target;removed_L+=l-new_l;n+=count
        e=target;l=new_l;t+=dt
        if gate_closes:status='EMISSION_GATE_CLOSED';break
    eligible_encounters=encounters
    if status in ('NO_RADIATION','EMISSION_GATE_CLOSED','CAPTURE_ANGULAR_THRESHOLD_NOT_MET'):
        extra=max(0.,horizon_ps-t)*2/clock.rates(e)['full_signed_period_ps']
        encounters+=extra
        if status!='EMISSION_GATE_CLOSED' and pericenter_gev_inverse(model.mass_gev,e,l,model.lambda_ev)<=gate:
            eligible_encounters+=extra
    survival=math.exp(-model.annihilation_optical_depth_per_eligible_encounter*eligible_encounters)
    excited_annihilation=1-survival
    capture_radiation=0.;capture_removed_L=0.;prob=excited_annihilation;pcap=0.;singlet_live=0.
    if capture_time is not None:
        capture_radiation=e+model.captured_binding_energy_gev
        capture_removed_L=l; e=-model.captured_binding_energy_gev;l=0.;pcap=survival
        age=max(0.,horizon_ps-capture_time)
        singlet_live=pcap*model.captured_singlet_fraction*math.exp(-age*model.singlet_annihilation_width_gev/HBAR_GEV_PS)
        prob+=pcap*model.captured_singlet_fraction-singlet_live
    energy_residual=initial_energy_gev-e-radiated-capture_radiation
    angular_residual=initial_L-l-removed_L-capture_removed_L
    if abs(energy_residual)>1e-9*max(1.,initial_energy_gev) or abs(angular_residual)>1e-9*max(1.,initial_L):
        raise RuntimeError('secular energy/angular ledger failed')
    other=pcap*(1-model.captured_singlet_fraction)
    uncaptured=survival if capture_time is None else 0.
    if abs(prob+other+singlet_live+uncaptured-1)>1e-12:
        raise RuntimeError('competing-outcome probability ledger failed')
    return {'scope':'DECLARED_SECULAR_CENTRAL_MODEL','status':status,'model':model.__dict__,
            'time_ps':t,'horizon_ps':horizon_ps,'capture_time_ps':capture_time,
            'initial_energy_gev':initial_energy_gev,'initial_L':initial_L,
            'final_energy_gev':e,'final_L':l,'mean_effective_emissions':n,
            'radial_encounters_in_secular_clock':encounters,
            'eligible_encounters_in_secular_clock':eligible_encounters,
            'excited_annihilation_probability_in_declared_model':excited_annihilation,
            'radiated_cooling_energy_gev':radiated,'capture_radiation_gev':capture_radiation,
            'orbital_L_removed_in_cooling':removed_L,'orbital_L_removed_at_capture':capture_removed_L,
            'energy_ledger_residual_gev':energy_residual,'angular_ledger_residual_hbar':angular_residual,
            'capture_probability_in_declared_model':pcap,
            'annihilation_probability_in_declared_model':prob if status!='SECULAR_UPDATE_OUTSIDE_ALLOWED_CENTRAL_ORBITS' else None,
            'probability_valid_through_ps':horizon_ps if status!='SECULAR_UPDATE_OUTSIDE_ALLOWED_CENTRAL_ORBITS' else t,
            'captured_other_spin_sector_probability':other,
            'captured_singlet_survival_probability':singlet_live,
            'uncaptured_survival_probability':uncaptured,
            'ledger_scope':'conditional surviving-action path; not survival-weighted ensemble radiation',
            'physical_prediction_status':'MODEL_DEPENDENT_NOT_NG_OR_QUANTUM_MATCHED',
            'history':history}
