"""Bare spin-resolved neutral scattering coefficients, without state occupations."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.neutral_hard_scattering import bare_neutral_kernels,spin_averages


def calculate():
    cases=[]
    for beta in (0.,.4):
        coarse=bare_neutral_kernels(225.,beta,2,.04,n_theta=24)
        fine=bare_neutral_kernels(225.,beta,2,.04,n_theta=48)
        channels={}
        for name,k in fine.items():
            error=float(np.max(abs(k-coarse[name]))/np.max(abs(k)))
            if error>1e-8:raise RuntimeError('angular quadrature did not converge')
            channels[name]={**spin_averages(k),'spin_matrix_real':k.real.tolist(),
                            'spin_matrix_imag':k.imag.tolist(),'relative_quadrature_change':error}
        cases.append({'beta_each_cm':beta,'v_moller':2*beta,'channels':channels})
    return {'scope':'BARE_ON_SHELL_NEUTRAL_DIRAC_TWO_BODY_SCATTERING',
            'input':{'mass_gev':225.,'n_ic':2,'alpha_ic_hard':.04,
                     'mass_convention':'degenerate tree-level neutral/charged fermion masses',
                     'spin_basis':'physical canonical particle/antiparticle spins, incoming momenta +/- z'},
            'cases':cases,'incoming_quantum_state':None,'annihilation_probability':None,
            'missing':['incoming spin and orbital density at the return','potential/IC-singlet wavefunction evolution',
                       'higher-order, higher-multiplicity and loop channels','radiation/capture dynamics']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=calculate();a.output.write_text(json.dumps(r,indent=2)+'\n')
    for c in r['cases']:
        print('beta',c['beta_each_cm'])
        for name,k in c['channels'].items():
            print(name,'singlet',k['singlet_sigma_v_gev_minus2'],'triplet-average',k['triplet_average_sigma_v_gev_minus2'])
