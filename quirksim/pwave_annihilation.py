"""Specified spin-triplet P-wave two-vector hard coefficients.

Cartesian orbital x spin basis; derivative of the full t/u amplitude at
threshold. Not an S-wave contact-density approximation or formation model.
"""
import math
import numpy as np
from scipy.sparse import bmat,diags
from scipy.sparse.linalg import eigsh
from .production_spin import GAMMA,u_spinors,v_spinors
from .neutral_hard_scattering import _slash
from .neutral_triplet import potential

TRIPLET=np.array([[-1,1j,0],[0,0,1],[0,0,1],[1,1j,0]],complex)/math.sqrt(2)


def p_projectors():
    swap=np.zeros((9,9));trace=np.eye(3).reshape(9)/math.sqrt(3)
    for i in range(3):
        for j in range(3):swap[3*i+j,3*j+i]=1
    p0=np.outer(trace,trace);p1=(np.eye(9)-swap)/2
    return {0:p0,1:p1,2:(np.eye(9)+swap)/2-p0}


def derivative_kernel(mass,vector_mass,alpha,identical=True,n_theta=12,n_phi=12,step_fraction=1e-4,vector_mass2=None,diagram='both'):
    mv2=vector_mass if vector_mass2 is None else vector_mass2
    if not all(math.isfinite(v) for v in (mass,vector_mass,mv2,alpha,step_fraction)):raise ValueError('finite inputs required')
    if 2*mass<=vector_mass+mv2 or mass<=0 or min(vector_mass,mv2)<0 or not 0<alpha<1 or not 0<step_fraction<.01 or diagram not in ('both','t','u'):
        raise ValueError('open channel and small positive threshold step required')
    if n_theta<4 or n_phi<8:raise ValueError('insufficient angular quadrature')
    step=mass*step_fraction;energy=math.hypot(mass,step);s=4*energy**2
    e1=(s+vector_mass**2-mv2**2)/(4*energy);e2=2*energy-e1;k=math.sqrt(e1**2-vector_mass**2)
    kernel=np.zeros((9,9),complex);nodes,weights=np.polynomial.legendre.leggauss(n_theta)
    for cosine,w in zip(nodes,weights):
        sine=math.sqrt(1-cosine*cosine)
        for phi in np.arange(n_phi)*2*math.pi/n_phi:
            direction=np.array([sine*math.cos(phi),sine*math.sin(phi),cosine])
            et=np.array([cosine*math.cos(phi),cosine*math.sin(phi),-sine]);ep=np.array([-math.sin(phi),math.cos(phi),0])
            k1=np.r_[e1,k*direction];k2=np.r_[e2,-k*direction]
            pol1=[np.r_[0,et],np.r_[0,ep]];pol2=pol1.copy()
            if vector_mass:
                pol1+=[np.r_[k/vector_mass,e1*direction/vector_mass]]
            if mv2:
                pol2+=[np.r_[k/mv2,-e2*direction/mv2]]
            for eps1 in pol1:
                for eps2 in pol2:
                    s1=_slash(eps1);s2=_slash(eps2);derivatives=[]
                    for axis in np.eye(3):
                        amplitudes=[]
                        for sign in (-1,1):
                            p=sign*step*axis;p1=np.r_[energy,p]
                            bar=v_spinors(-p,mass).conj().T@GAMMA[0]
                            propagators=[]
                            for kk in (k1,k2):
                                q=p1-kk;propagators.append((_slash(q)+mass*np.eye(4))/(q[0]**2-q[1:]@q[1:]-mass**2))
                            op=4*math.pi*alpha*((s2@propagators[0]@s1 if diagram!='u' else np.zeros((4,4)))+
                                (s1@propagators[1]@s2 if diagram!='t' else np.zeros((4,4))))
                            amplitude=(bar@op@u_spinors(p,mass)).T.reshape(4)
                            amplitudes.append(TRIPLET.T@amplitude)
                        derivatives.append((amplitudes[1]-amplitudes[0])/(2*step))
                    a=np.array(derivatives).reshape(9)
                    kernel+=w/(2*n_phi)*np.outer(a.conj(),a)
    return kernel*(k/energy)/(32*math.pi*energy**2)/(2 if identical else 1)


def j_averages(kernel):
    return {str(j):float(np.trace(p@kernel).real/(2*j+1)) for j,p in p_projectors().items()}


def compact_p_state(mass,delta,n_ic,alpha,lambda_ev,points=2000):
    if (type(points) is not int or points<100 or type(n_ic) is not int or n_ic<2
        or not all(math.isfinite(v) for v in (mass,delta,alpha,lambda_ev))
        or mass<=0 or delta<0 or not 0<alpha<1 or lambda_ev<=0):raise ValueError('physical radial inputs required')
    cf=(n_ic*n_ic-1)/(2*n_ic);kappa=cf*alpha
    rmax=40/(mass*kappa/2);dr=rmax/(points+1);r=dr*np.arange(1,points+1)
    v=potential(r,delta,n_ic,alpha,lambda_ev)
    kinetic=diags([-np.ones(points-1),2*np.ones(points),-np.ones(points-1)],[-1,0,1])/(mass*dr**2)+diags(2/(mass*r*r))
    blocks=[[diags(v[:,i,j])+(kinetic if i==j else 0) for j in range(3)] for i in range(3)]
    h=bmat(blocks,format='csc');energy,u=eigsh(h,k=1,sigma=-mass*(kappa+.1)**2,which='LM',tol=1e-11,v0=np.ones(3*points))
    u=u[:,0].reshape(3,points)/math.sqrt(dr)
    if u[1,0]<0:u=-u
    derivative=[np.polynomial.polynomial.polyfit(r[:5],row[:5]/r[:5]**2,2)[0] for row in u]
    return {'energy_gev':float(energy[0]),'radial_derivative_at_origin_gev_5halves':derivative,
        'points':points,'dr_gev_inverse':dr,'flavor_probabilities':(np.sum(u*u,axis=1)*dr).tolist()}
