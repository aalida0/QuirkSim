"""Specified LO Dirac Y=0 triplet, IC-singlet, spin-zero quantum sector.

Basis: (Q+ anti-Q+, Q0 anti-Q0, Q- anti-Q-). This is not the two-channel
Majorana-wino basis. Amplitude phases put W-exchange off-diagonal potentials
negative. Formation and the spin-one sector are not specified here.
"""
import math
import numpy as np
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import eigsh
from .physics import GF, MW, MZ


def ew_inputs():
    aw=math.sqrt(2)*GF*MW**2/math.pi
    sw2=1-(MW/MZ)**2
    return aw, sw2, aw*sw2


def hard_spin_zero(mass_gev, n_ic, alpha_ic_hard):
    """Hermitian sigma*v matrices; width = psi(0)^dagger K psi(0).

    Strict threshold and degenerate tree-level mediator masses. The radiative
    splitting is used in the potential thresholds, not inserted into only
    selected propagators of a tree-level gauge amplitude.
    """
    if (not math.isfinite(mass_gev) or mass_gev<=0 or type(n_ic) is not int
            or n_ic<2 or not math.isfinite(alpha_ic_hard) or not 0<=alpha_ic_hard<1):
        raise ValueError('positive mass, integer N>=2, perturbative coupling required')
    aw,s,ae=ew_inputs();c=1-s;m=mass_gev
    def phase(mv):
        r=(mv/m)**2
        return (1-r)**1.5/(1-r/2)**2 if r<1 else 0.
    def matrix(coefficient,amplitude):
        v=np.array(amplitude,float);return coefficient*np.outer(v,v)
    base=math.pi/m**2
    return {
        'ic_gluon_pair':matrix(base*alpha_ic_hard**2*(n_ic*n_ic-1)/n_ic,[1,1,1]),
        'W_plus_W_minus':matrix(8*base*n_ic*aw**2*phase(MW),[.5,1,.5]),
        'gamma_gamma':matrix(4*base*n_ic*ae**2,[1,0,1]),
        'gamma_Z':matrix(8*base*n_ic*ae*aw*c*max(0.,1-(MZ/(2*m))**2),[1,0,1]),
        'Z_Z':matrix(4*base*n_ic*(aw*c)**2*phase(MZ),[1,0,1])}


def potential(r, delta_gev, n_ic, alpha_ic_potential, lambda_ev):
    """Static LO central potential in GeV; r in GeV^-1.

    Single gauge exchange derived from triplet generators, plus declared
    fixed-coupling Cornell IC potential. No nonperturbative emission law.
    """
    r=np.asarray(r,float)
    if (r.ndim!=1 or not np.isfinite(r).all() or np.any(r<=0)
            or not all(map(math.isfinite,[delta_gev,alpha_ic_potential,lambda_ev]))
            or delta_gev<0 or not 0<alpha_ic_potential<1 or lambda_ev<0
            or type(n_ic) is not int or n_ic<2):
        raise ValueError('physical potential parameters and positive radii required')
    aw,s,ae=ew_inputs();cf=(n_ic*n_ic-1)/(2*n_ic)
    common=-cf*alpha_ic_potential/r+(lambda_ev*1e-9)**2*r
    a=(ae+aw*(1-s)*np.exp(-MZ*r))/r
    b=aw*np.exp(-MW*r)/r
    v=np.zeros((len(r),3,3))
    v[:,0,0]=common+2*delta_gev-a;v[:,1,1]=common;v[:,2,2]=v[:,0,0]
    v[:,0,1]=v[:,1,0]=-b;v[:,1,2]=v[:,2,1]=-b
    return v


def ground_state(mass_gev, delta_gev, n_ic, alpha_ic_potential, lambda_ev, points=4000):
    """Normalized spin-independent l=0 radial eigenstate in a fixed model.

    u(0)=u(rmax)=0, integral dr sum |u_i|^2=1. Check mesh convergence before
    using contact amplitudes. This does not infer occupation from a trajectory.
    """
    if (not math.isfinite(mass_gev) or mass_gev<=0 or type(points) is not int or points<100
            or type(n_ic) is not int or n_ic<2 or not math.isfinite(alpha_ic_potential)
            or not 0<alpha_ic_potential<1):
        raise ValueError('positive mass, perturbative coupling, N>=2 and >=100 mesh points required')
    cf=(n_ic*n_ic-1)/(2*n_ic);kappa=cf*alpha_ic_potential
    if kappa<=0:raise ValueError('positive Coulomb attraction required')
    rmax=30/(mass_gev*kappa/2);dr=rmax/(points+1);r=dr*np.arange(1,points+1)
    v=potential(r,delta_gev,n_ic,alpha_ic_potential,lambda_ev)
    kinetic=diags([-np.ones(points-1),2*np.ones(points),-np.ones(points-1)],[-1,0,1])/(mass_gev*dr*dr)
    blocks=[[diags(v[:,i,j])+(kinetic if i==j else 0) for j in range(3)] for i in range(3)]
    h=bmat(blocks,format='csc');aw,_,ae=ew_inputs()
    shift=-mass_gev*(kappa+2*aw+ae)**2
    # Shift is below the Coulomb lower bound; nearest eigenvalue is the ground.
    energy,u=eigsh(h,k=1,sigma=shift,which='LM',tol=1e-11,v0=np.ones(3*points))
    u=u[:,0].reshape(3,points)/math.sqrt(dr)
    if u[1,0]<0:u=-u
    psi=np.array([np.polynomial.polynomial.polyfit(r[:5],row[:5]/r[:5],2)[0]
                  for row in u])/math.sqrt(4*math.pi)
    return {'energy_relative_to_neutral_threshold_gev':float(energy[0]),
            'psi_at_origin_gev_three_halves':psi.tolist(),
            'flavor_probabilities':(np.sum(u*u,axis=1)*dr).tolist(),
            'points':points,'dr_gev_minus1':dr,'rmax_gev_minus1':rmax,
            'occupation_in_transport':None}


def widths_from_state(contact_amplitudes, kernels):
    psi=np.asarray(contact_amplitudes,complex)
    if psi.shape!=(3,) or not np.isfinite(psi).all():raise ValueError('three finite contact amplitudes required')
    result={}
    for name,k in kernels.items():
        k=np.asarray(k,complex)
        if (k.shape!=(3,3) or not np.isfinite(k).all() or np.max(abs(k-k.conj().T))>1e-12
                or np.linalg.eigvalsh(k).min() < -1e-12):raise ValueError('positive Hermitian kernel required')
        result[name]=float(np.real(np.vdot(psi,k@psi)))
    return result


def gaussian_overlap_log_ratio(impact_um, transverse_sigma_um):
    """Diagnostic ratio in a supplied ballistic Gaussian contact ansatz.

    Width is not fixed by mean classical coordinates. Ignores capture,
    radiation, focusing and spin. This is a non-identifiability demonstration,
    not the triplet physical suppression factor.
    """
    if not (math.isfinite(impact_um) and impact_um>=0 and math.isfinite(transverse_sigma_um) and transverse_sigma_um>0):
        raise ValueError('nonnegative impact and positive packet width required')
    return -.5*(impact_um/transverse_sigma_um)**2
