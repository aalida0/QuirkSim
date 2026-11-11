"""Reproducible compact-state partial annihilation benchmark, independent of formation."""
import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.physics import benchmark
from quirksim.radiation_annihilation import (
    ALPHA_ZERO, coulomb_ground_reference, spin_singlet_two_boson_widths)


def calculate():
    m=benchmark('triplet',225.)['charged_mass_gev']
    n=2; alpha_potential=.04; alpha_hard=.04
    cf=(n*n-1)/(2*n)
    # Explicit single-channel fixed-coupling approximation, including photon.
    kappa=cf*alpha_potential+ALPHA_ZERO
    reference=coulomb_ground_reference(m,m,kappa)
    widths=spin_singlet_two_boson_widths(m,reference['psi_squared_gev3'],n,alpha_hard)
    mu=m/2; radius=1/(mu*kappa); tension=(774e-9)**2
    return {'scope':'COMPACT_CHARGED_COMPONENT_LO_PARTIAL_WIDTH_BENCHMARK',
            'triplet_transport_annihilation_probability':None,
            'input':{'neutral_mass_gev':225.,'charged_mass_gev':m,'n_ic':n,
                     'alpha_ic_potential':alpha_potential,'alpha_ic_hard':alpha_hard,
                     'alpha_em':ALPHA_ZERO,'lambda_tension_ev':774.,
                     'spin':0,'orbital_l':0,'radial_n':1,
                     'potential':'-(C_F alpha_ic_potential + alpha_em)/r',
                     'coupling_scheme':'fixed-coupling LO illustration; no Lambda_MS matching'},
            'state':reference,'partial_annihilation':widths,
            'scales':{'bohr_momentum_gev':1/radius,
                      'mean_linear_potential_gev':1.5*tension*radius,
                      'linear_to_binding_ratio':1.5*tension*radius/reference['binding_energy_gev']},
            'omissions':['finite-range Z potential and W-induced flavor mixing',
                         'running of the potential coupling',
                         'radiative capture and weak-decay competition in the compact state'],
            'release':'PARTIAL_WIDTH_FORMULA_IMPLEMENTED; SINGLE_CHANNEL_NUMBERS_APPROXIMATE; TRANSPORT_MATCH_UNKNOWN'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=calculate();a.output.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r['partial_annihilation'],indent=2))
