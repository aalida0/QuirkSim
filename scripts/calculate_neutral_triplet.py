"""Specified neutral-connected spin-zero eigenstate and matching ambiguity audit."""
import argparse
import json
import math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.physics import benchmark, HBAR_GEV_PS
from quirksim.neutral_triplet import ground_state, hard_spin_zero, widths_from_state, gaussian_overlap_log_ratio, ew_inputs


def calculate():
    mass=225.;delta=benchmark('triplet',mass)['charged_mass_gev']-mass
    kernels=hard_spin_zero(mass,2,.04);meshes=[]
    for points in (4000,8000,16000,32000):
        state=ground_state(mass,delta,2,.04,774.,points=points)
        widths=widths_from_state(state['psi_at_origin_gev_three_halves'],kernels)
        meshes.append({'state':state,'partial_widths_gev':widths})
    final=meshes[-1];previous=meshes[-2]
    relative={k:abs(v/previous['partial_widths_gev'][k]-1) for k,v in final['partial_widths_gev'].items()}
    if max(relative.values())>1e-3:raise RuntimeError('contact widths not converged to 0.1 percent between final meshes')
    total=sum(final['partial_widths_gev'].values())
    impact=4.154394665804169
    return {'scope':'SPECIFIED_LO_SPIN_ZERO_COMPACT_STATE_NOT_TRANSPORT_ANNIHILATION',
            'model':{'neutral_mass_gev':mass,'charged_minus_neutral_gev':delta,
                     'n_ic':2,'alpha_ic_potential':.04,'alpha_ic_hard':.04,'lambda_tension_ev':774.,
                     'representation':'DIRAC_SU2_TRIPLET_Y0_IC_FUNDAMENTAL',
                     'ew_scheme':'tree GF MW MZ; alpha_EM derived consistently',
                     'ew_alphaW_sin2W_alphaEM':ew_inputs(),
                     'potential':'fixed-coupling Cornell IC plus single-exchange EW matrix',
                     'flavor_basis':['Q+ antiQ+','Q0 antiQ0','Q- antiQ-']},
            'mesh_checks':meshes,'last_relative_width_changes':relative,
            'sum_modeled_two_boson_widths_gev':total,'hbar_over_modeled_sum_ps':HBAR_GEV_PS/total,
            'transport_matching':{'status':'NOT_IDENTIFIABLE_FROM_CLASSICAL_SAMPLE_MEANS',
                                  'spin_density_matrix':None,'wavepacket_density_operator':None,
                                  'formation_probability':None,'radiation_transition_kernel':None},
            'gaussian_counterexample':{
                'scope':'BALLISTIC_GAUSSIAN_CONTACT_ANSATZ_ONLY_NOT_A_PHYSICAL_ERROR_BAND',
                'recoil_impact_um':impact,'reference_impact_um':0.,
                'same_means_different_widths':[
                    {'transverse_sigma_um':s,'log_recoil_to_reference_contact_ratio':gaussian_overlap_log_ratio(impact,s)}
                    for s in (.000255,.001,.1,1.,10.)]},
            'physical_annihilation_probability_with_recoil':None,
            'physical_annihilation_probability_without_recoil':None,
            'missing':['spin-one annihilation and transitions','formation from production and weak decay',
                       'finite-energy and higher-order hard corrections','nonperturbative radiation and recoil',
                       'hard/Bohr running-coupling and string-tension scale matching'],
            'release':'COMPACT_SPIN_ZERO_MODEL_CALCULATION_ONLY; ITEMS_2_AND_3_UNRESOLVED'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=calculate();a.output.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({'widths_kev':{k:v*1e6 for k,v in r['mesh_checks'][-1]['partial_widths_gev'].items()},
                      'last_relative_width_changes':r['last_relative_width_changes'],
                      'hbar_over_modeled_sum_ps':r['hbar_over_modeled_sum_ps']},indent=2))
