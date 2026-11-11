"""Canonical Dirac spins for a specified vector-current production model.

These local amplitudes do not supply continuous spin transport in a string.
Basis is (up up, up down, down up, down down), particle then antiparticle.
"""
import math
import numpy as np

PAULI=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
EPS=np.array([[0,-1],[1,0]],complex)
ZERO=np.zeros((2,2),complex); ID=np.eye(2,dtype=complex)
GAMMA=[np.block([[ID,ZERO],[ZERO,-ID]])]+[np.block([[ZERO,s],[-s,ZERO]]) for s in PAULI]
SINGLET=np.array([0,1,-1,0],complex)/math.sqrt(2)


def u_spinors(momentum,mass):
    p=np.asarray(momentum,float)
    if p.shape!=(3,) or not np.isfinite(p).all() or not math.isfinite(mass) or mass<=0:
        raise ValueError('finite momentum and positive mass required')
    e=math.hypot(mass,float(np.linalg.norm(p)));a=math.sqrt(e+mass)
    return np.vstack([a*ID,np.einsum('i,ijk->jk',p,PAULI)/a])


def v_spinors(momentum,mass):
    u=u_spinors(momentum,mass)
    return np.vstack([u[2:]@EPS,u[:2]@EPS])


def vector_currents(momentum,mass,beam=(0.,0.,1.)):
    """Two circular transverse source amplitudes for massless incoming q qbar.

    The final heavy current is vectorlike (gamma/Z charged-triplet source).
    Initial helicity weights encode quark couplings and luminosities; they
    do not affect the absence of a CM singlet in this production amplitude.
    """
    beam=np.asarray(beam,float)
    if beam.shape!=(3,) or not np.isfinite(beam).all() or np.linalg.norm(beam)==0:
        raise ValueError('finite nonzero beam direction required')
    beam=beam/np.linalg.norm(beam);axis=np.eye(3)[int(np.argmin(abs(beam)))]
    ex=np.cross(beam,axis);ex/=np.linalg.norm(ex);ey=np.cross(beam,ex)
    u=u_spinors(momentum,mass);v=v_spinors(-np.asarray(momentum,float),mass)
    current=np.array([(u.conj().T@GAMMA[0]@g@v).reshape(4) for g in GAMMA[1:]])
    jx=ex@current;jy=ey@current
    return np.array([(jx+1j*jy)/math.sqrt(2),(jx-1j*jy)/math.sqrt(2)])


def dy_vector_density(momentum,mass,beam=(0.,0.,1.),helicity_weights=(.5,.5)):
    weights=np.asarray(helicity_weights,float)
    if weights.shape!=(2,) or not np.isfinite(weights).all() or np.any(weights<0) or weights.sum()<=0:
        raise ValueError('nonnegative nonzero helicity weights required')
    amplitudes=vector_currents(momentum,mass,beam)
    rho=sum(w*np.outer(a,a.conj()) for w,a in zip(weights,amplitudes))
    return rho/np.trace(rho).real


def singlet_fraction(rho):
    rho=np.asarray(rho,complex)
    if (rho.shape!=(4,4) or not np.isfinite(rho).all() or np.max(abs(rho-rho.conj().T))>1e-12
            or abs(np.trace(rho)-1)>1e-12 or np.linalg.eigvalsh(rho).min() < -1e-12):
        raise ValueError('normalized positive spin density required')
    return float(np.vdot(SINGLET,rho@SINGLET).real)


def pion_spin_map(parent_momentum,parent_mass,daughter_momentum,daughter_mass):
    """Normalized canonical-spin map of the LO vector-current pion vertex.

    q_mu ubar_d gamma^mu u_P=(M-m) ubar_d u_P. At fixed momenta the 2x2
    scalar-bilinear map has M^dagger M proportional to identity. It preserves
    spin exactly in the parent-rest canonical basis; other frames induce a
    canonical Wigner rotation. The physical antiparticle map is the same up
    to an irrelevant common phase. No propagation between vertices is implied.
    """
    if not parent_mass>daughter_mass:raise ValueError('ordered masses required')
    parent=u_spinors(parent_momentum,parent_mass);daughter=u_spinors(daughter_momentum,daughter_mass)
    amplitude=daughter.conj().T@GAMMA[0]@parent
    scale=float(np.trace(amplitude.conj().T@amplitude).real/2)
    rotation=amplitude/math.sqrt(scale)
    if np.max(abs(rotation.conj().T@rotation-ID))>1e-10:raise RuntimeError('pion map not unitary')
    return rotation


def canonical_boost_rotation(momentum,mass,beta):
    """Wigner rotation from a rotation-free active Lorentz boost."""
    from .decays import boost
    b=np.asarray(beta,float)
    if b.shape!=(3,) or not np.isfinite(b).all() or b@b>=1:raise ValueError('finite subluminal boost required')
    if b@b==0:return ID.copy()
    gamma=1/math.sqrt(1-float(b@b));n=b/np.linalg.norm(b)
    alpha=sum(n[j]*(GAMMA[0]@GAMMA[j+1]) for j in range(3))
    lorentz=math.sqrt((gamma+1)/2)*np.eye(4)+math.sqrt((gamma-1)/2)*alpha
    source=u_spinors(momentum,mass)
    target_p=boost(math.hypot(mass,np.linalg.norm(momentum)),momentum,b)[1:]
    target=u_spinors(target_p,mass)
    w=target.conj().T@GAMMA[0]@lorentz@source/(2*mass)
    if np.max(abs(w.conj().T@w-ID))>1e-10:raise RuntimeError('boost rotation not unitary')
    return w


def rotate_pair(rho,rotation0,rotation1):
    singlet_fraction(rho)  # validate the density, independent of its fraction
    for w in (rotation0,rotation1):
        if np.shape(w)!=(2,2) or not np.isfinite(w).all() or np.max(abs(w.conj().T@w-ID))>1e-10:
            raise ValueError('unitary local spin maps required')
    w=np.kron(rotation0,rotation1)
    return w@rho@w.conj().T
