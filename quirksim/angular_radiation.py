"""Necessary low-L reachability under declared finite-energy emission laws.

Kang/Luty motivates O(Lambda) energy and O(1) angular steps, not the hard
bounds below. No emitted spectrum, rate, capture or annihilation is inferred.
"""
import math


def reachability(initial_L,available_energy_gev,min_emission_energy_gev,
                 max_angular_step=1.,capture_max_L=1.):
    values=(initial_L,available_energy_gev,min_emission_energy_gev,max_angular_step,capture_max_L)
    if not all(math.isfinite(x) for x in values) or min(initial_L,available_energy_gev,max_angular_step,capture_max_L)<0 or min_emission_energy_gev<=0:
        raise ValueError('finite nonnegative actions/energy and positive energy floor required')
    n=math.floor(available_energy_gev/min_emission_energy_gev)
    distance=max(0.,initial_L-capture_max_L)
    possible=distance<=n*max_angular_step
    # ||sum Delta J||^2 is a nonnegative submartingale for zero conditional
    # vector drift. Its expectation is <= N*jmax^2. The maximal inequality
    # bounds hitting the low-L ball at ANY emission, not just at the last.
    upper=(1. if distance==0 else min(1.,n*(max_angular_step/distance)**2)) if possible else 0.
    return {'max_emissions':n,'max_directed_angular_removal':n*max_angular_step,
        'angular_gap_to_capture':distance,'reachable_by_any_bounded_step_sequence':possible,
        'unbiased_vector_walk_low_L_hitting_probability_upper_bound':upper,
        'deterministic_remaining_L_lower_bound':max(0.,initial_L-n*max_angular_step),
        'scope':'CONDITIONAL_NECESSARY_LOW_L_REACHABILITY_NOT_ANNIHILATION'}
