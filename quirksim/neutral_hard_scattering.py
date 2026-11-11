"""Finite-velocity bare neutral Dirac-pair two-vector annihilation kernels.

No potential resummation, incoming wavepacket, spin transport or formation.
Incoming momenta are +/-z in their CM frame. Sigma*v uses v_Moller=2 beta.
"""
import math
import numpy as np
from .production_spin import GAMMA,u_spinors,v_spinors,SINGLET
from .neutral_triplet import ew_inputs
from .physics import MW


def _slash(p):
    return GAMMA[0]*p[0]-sum((GAMMA[j+1]*p[j+1] for j in range(3)),np.zeros((4,4),complex))


def equal_vector_kernel(mass,beta,vector_mass,alpha,identical,n_theta=32,n_phi=12):
    """Vectorlike t/u amplitudes with degenerate mediator and external masses.

    Applicable to neutral-triplet WW (the neutral s-channel is absent) and
    the abelian kinematic part of IC-singlet two-gluon annihilation.
    Returns a 4x4 physical canonical-spin kernel in GeV^-2.
    """
    if (not all(map(math.isfinite,[mass,beta,vector_mass,alpha])) or mass<=0
            or not 0<=beta<1 or vector_mass<0 or not 0<alpha<1
            or type(n_theta) is not int or n_theta<4 or type(n_phi) is not int or n_phi<8):
        raise ValueError('physical masses/coupling, beta<1 and sufficient quadrature required')
    energy=mass/math.sqrt(1-beta*beta);p=beta*energy
    if energy<=vector_mass:return np.zeros((4,4),complex)
    momentum=math.sqrt(energy*energy-vector_mass*vector_mass)
    p1=np.array([energy,0,0,p]);particle=u_spinors(p1[1:],mass)
    anti=v_spinors(-p1[1:],mass);bar=anti.conj().T@GAMMA[0]
    nodes,weights=np.polynomial.legendre.leggauss(n_theta)
    kernel=np.zeros((4,4),complex)
    for cosine,weight in zip(nodes,weights):
        sine=math.sqrt(max(0.,1-cosine*cosine))
        for phi in np.arange(n_phi)*2*math.pi/n_phi:
            n=np.array([sine*math.cos(phi),sine*math.sin(phi),cosine])
            et=np.array([cosine*math.cos(phi),cosine*math.sin(phi),-sine])
            ep=np.array([-math.sin(phi),math.cos(phi),0])
            k1=np.r_[energy,momentum*n];k2=np.r_[energy,-momentum*n]
            polar1=[np.r_[0,et],np.r_[0,ep]];polar2=polar1.copy()
            if vector_mass:
                polar1+=[np.r_[momentum/vector_mass,energy*n/vector_mass]]
                polar2+=[np.r_[momentum/vector_mass,-energy*n/vector_mass]]
            q1=p1-k1;q2=p1-k2
            d1=q1[0]**2-q1[1:]@q1[1:]-mass*mass
            d2=q2[0]**2-q2[1:]@q2[1:]-mass*mass
            prop1=(_slash(q1)+mass*np.eye(4))/d1;prop2=(_slash(q2)+mass*np.eye(4))/d2
            for eps1 in polar1:
                s1=_slash(eps1)
                for eps2 in polar2:
                    s2=_slash(eps2)
                    op=(s2@prop1@s1+s1@prop2@s2)*(4*math.pi*alpha)
                    # M[p-spin, anti-spin]; the vbar matrix has anti-spin rows.
                    amplitude=(bar@op@particle).T.reshape(4)
                    kernel+=weight/(2*n_phi)*np.outer(amplitude.conj(),amplitude)
    kernel*=momentum/energy/(32*math.pi*energy*energy)/(2 if identical else 1)
    return kernel


def bare_neutral_kernels(mass,beta,n_ic,alpha_ic,n_theta=32):
    if type(n_ic) is not int or n_ic<2 or not math.isfinite(alpha_ic) or not 0<alpha_ic<1:
        raise ValueError('N>=2 and perturbative IC coupling required')
    aw,_,_=ew_inputs()
    # EW amplitudes have a sqrt(N) IC-singlet projection. For two IC gluons,
    # the ratio to abelian same-alpha photons is (N^2-1)/(4N).
    return {'W_plus_W_minus':n_ic*equal_vector_kernel(mass,beta,MW,aw,False,n_theta),
            'ic_gluon_pair':(n_ic*n_ic-1)/(4*n_ic)*equal_vector_kernel(mass,beta,0.,alpha_ic,True,n_theta)}


def spin_averages(kernel):
    singlet=float(np.real(np.vdot(SINGLET,kernel@SINGLET)))
    return {'singlet_sigma_v_gev_minus2':singlet,
            'triplet_average_sigma_v_gev_minus2':float((np.trace(kernel).real-singlet)/3),
            'unpolarized_sigma_v_gev_minus2':float(np.trace(kernel).real/4)}
