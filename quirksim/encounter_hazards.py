"""Uncertainty-bracketed *diagnostic* encounter, radiation, annihilation hazards.

NOT physical quirk annihilation or a self-consistent radiation dynamics model.
Published phenom often assumes a glueball emitted with probability epsilon
per oscillation, of energy ~Lambda_IC (Li et al., PRD 109, 095005, 2024).
Its absolute rate is unknown. Quantum annihilation cannot be reduced to a
classical distance threshold.

This module intentionally DOES NOT alter the string trajectory or pair energy.
It evaluates multiple explicitly user-selected independent-hazard scenarios
and exposes pre-decay survivorship biases. Any result is conditional.
"""
from dataclasses import dataclass,replace
import math

from .dynamics import HBARC_GEV_CM


@dataclass(frozen=True)
class HazardModel:
    lambda_ev: float
    quirk_mass_gev: float
    annihilation_per_geometric_contact: float = 0.0
    glueball_emission_per_contact: float = 0.0
    glueball_energy_in_units_lambda: float = 1.0
    radius_definition: str = 'one_over_lambda'

    def __post_init__(self):
        if self.lambda_ev<=0 or self.quirk_mass_gev<=0:
            raise ValueError('positive Lambda and mass required')
        if not 0<=self.annihilation_per_geometric_contact<=1:
            raise ValueError('invalid annihilation hazard')
        if not 0<=self.glueball_emission_per_contact<=1:
            raise ValueError('invalid emission probability')
        if self.glueball_energy_in_units_lambda<=0:
            raise ValueError('invalid glueball energy scale')
        if self.radius_definition not in ('one_over_lambda','one_over_mass','geometric_mean'):
            raise ValueError('unknown reference radius')

    @property
    def radius_cm(self):
        lam=self.lambda_ev*1e-9
        if self.radius_definition=='one_over_lambda':
            ref=lam
        elif self.radius_definition=='one_over_mass':
            ref=self.quirk_mass_gev
        else:
            ref=math.sqrt(self.quirk_mass_gev*lam)
        return HBARC_GEV_CM/ref


def evaluate(encounters, first_weak_decay_ps, model):
    """Evaluate a conditional scenario without integrating energy loss.

    Geometric encounters are a sequence of event-local minima with
    fields time_ps and separation_cm. Surviving pre-weak-decay contacts
    is explicitly required for the later weak decay to remain observable.
    """
    ordered=sorted(encounters,key=lambda x:x['time_ps'])
    if any(e['separation_cm']<0 or not math.isfinite(e['time_ps'])
           or not math.isfinite(e['separation_cm']) for e in ordered):
        raise ValueError('invalid encounter record')
    survival=1.0
    survival_to_weak=1.0
    n_contact=0
    n_pre=0
    expected_radiation=0.0
    weak_reached=False
    for e in ordered:
        t=e['time_ps']
        if first_weak_decay_ps is not None and not weak_reached and t>=first_weak_decay_ps:
            survival_to_weak=survival
            weak_reached=True
        if e['separation_cm']>model.radius_cm:
            continue
        n_contact+=1
        if first_weak_decay_ps is None or t<first_weak_decay_ps:
            n_pre+=1
        # Illustrative independent Bernoulli model: not a QCD prediction.
        expected_radiation+=survival*model.glueball_emission_per_contact*(
            model.glueball_energy_in_units_lambda*model.lambda_ev*1e-9)
        survival*=1-model.annihilation_per_geometric_contact
    if first_weak_decay_ps is not None and not weak_reached:
        survival_to_weak=survival
    return {
        'radius_definition':model.radius_definition,
        'radius_cm':model.radius_cm,
        'contacts':n_contact,'pre_weak_decay_contacts':n_pre,
        'conditional_survival_to_first_weak_decay':(
            survival_to_weak if first_weak_decay_ps is not None else None),
        'conditional_survival_to_observation_end':survival,
        'expected_glueball_radiation_gev_without_backreaction':expected_radiation,
        'p_ann_per_contact_assumed':model.annihilation_per_geometric_contact,
        'p_glue_per_contact_assumed':model.glueball_emission_per_contact,
        'status':'ASSUMPTION_ENVELOPE_ONLY',
        'warnings':[
           'geometric minima are not measured annihilation probabilities',
           'annihilation requires quantum short-distance wavefunctions',
           'glueball probability and energy are nonperturbative assumptions',
           'radiation energy and momentum NOT fed back into string state',
           'independence of successive contacts is an unvalidated assumption']
    }


def scenario_envelope(encounters,first_weak_decay_ps,lambda_ev,mass_gev):
    """Display deliberate extreme assumptions to expose non-identifiability."""
    result=[]
    for radius in ('one_over_mass','geometric_mean','one_over_lambda'):
        for p_ann in (0.,.01,.1,1.):
            for eps in (0.,.5,1.):
                model=HazardModel(lambda_ev,mass_gev,p_ann,eps,
                                  radius_definition=radius)
                result.append(evaluate(encounters,first_weak_decay_ps,model))
    return result
