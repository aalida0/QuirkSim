"""LO widths of specified lowest neutral-connected 3P_J compact states."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.pwave_annihilation import compact_p_state,derivative_kernel,j_averages
from quirksim.quirk_g_parity import invariant_checks,neutral_capture_channels
from quirksim.neutral_triplet import ew_inputs
from quirksim.physics import MW,MZ,HBAR_GEV_PS


def calculate(compact,grids=(4000,8000,16000)):
    data=json.loads(compact.read_text());m=data['model'];mass=m['neutral_mass_gev'];n=m['n_ic'];aw,s,ae=ew_inputs()
    if n!=2 or m['representation']!='DIRAC_SU2_TRIPLET_Y0_IC_FUNDAMENTAL':raise ValueError('G projection implemented only for the specified two-color Dirac triplet')
    specs={'ic_gluon_pair':(0.,0.,m['alpha_ic_hard'],True,(n*n-1)/(4*n),[1,1,1]),
        'W_plus_W_minus':(MW,MW,aw,False,n,[.5,1,.5]),
        'gamma_gamma':(0.,0.,ae,True,n,[1,0,1]),
        'gamma_Z':(0.,MZ,math.sqrt(ae*aw*(1-s)),False,n,[1,0,1]),
        'Z_Z':(MZ,MZ,aw*(1-s),True,n,[1,0,1])}
    kernels={};checks=[]
    for name,(v1,v2,alpha,identical,color,flavor) in specs.items():
        coarse=derivative_kernel(mass,v1,alpha,identical,12,12,1e-4,v2)
        fine=derivative_kernel(mass,v1,alpha,identical,20,16,5e-5,v2)
        scale=max(np.max(abs(fine)),1e-100)
        difference=float(np.max(abs(coarse-fine))/scale)
        if difference>2e-6:raise ValueError('hard coefficient quadrature/step convergence failed')
        if np.linalg.eigvalsh(fine).min()<-1e-9*scale:raise ValueError('nonpositive P-wave kernel')
        kernels[name]=(fine,color,flavor)
        checks.append({'channel':name,'relative_matrix_change':difference,'J_averaged_derivative_kernel_gev_minus4':j_averages(fine)})
    states=[]
    for points in grids:
        state=compact_p_state(mass,m['charged_minus_neutral_gev'],n,m['alpha_ic_potential'],m['lambda_tension_ev'],points)
        radial=np.array(state['radial_derivative_at_origin_gev_5halves']);widths={str(j):{} for j in (0,1,2)}
        for name,(kernel,color,flavor) in kernels.items():
            factor=3/(4*math.pi)*color*float(radial@flavor)**2
            for j,k in j_averages(kernel).items():
                widths[j][name]=0. if j=='1' and specs[name][0]==specs[name][1]==0. else max(0.,factor*k)
        state['partial_widths_gev']=widths
        state['sum_modeled_widths_gev']={j:sum(w.values()) for j,w in widths.items()}
        state['proper_lifetime_from_modeled_channels_ps']={j:HBAR_GEV_PS/w for j,w in state['sum_modeled_widths_gev'].items()}
        states.append(state)
    shifts={j:abs(states[-1]['sum_modeled_widths_gev'][j]/states[-2]['sum_modeled_widths_gev'][j]-1) for j in ('0','1','2')}
    if max(shifts.values())>.01:raise ValueError('P-state width mesh change exceeds one percent')
    return {'scope':'SPECIFIED_LO_G_EVEN_COMPACT_P_STATES_NOT_FORMATION_PROBABILITY',
        'compact_input_sha256':hashlib.sha256(compact.read_bytes()).hexdigest(),'model':m,
        'twisted_charge_conjugation_checks':invariant_checks(),'neutral_source_G_even_selection':neutral_capture_channels(),
        'hard_coefficient_checks':checks,'states':states,'last_mesh_relative_width_change':shifts,
        'limitations':['Derivative-matched LO two-vector channels; higher multiplicity and radiative corrections omitted.',
            'Degenerate hard masses and fixed-coupling spin-independent potential; splitting only in potential thresholds.',
            'P-state spin-orbit/hyperfine splittings and capture occupations not calculated.',
            'G-parity requires minimal Y=0 vectorlike couplings and no G-breaking operators.',
            'The AF gun is a classical branch, not a coherent source amplitude; G selection is a specified physical production model.']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--compact',type=Path,default=Path(__file__).resolve().parents[1]/'studies/data/neutral_triplet_quantum_2026-11-11.json')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=calculate(a.compact)
    a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print('P-wave modeled widths GeV',r['states'][-1]['sum_modeled_widths_gev'])
    print('lifetimes ps',r['states'][-1]['proper_lifetime_from_modeled_channels_ps'])
