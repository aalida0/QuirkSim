"""Independent higher-order spatial quadrature for discrete Nambu--Goto action.

A finite element has piecewise-linear x(u) and v(u), but the Nambu--Goto
velocity factor sqrt(1-v_perp(u)^2) is NOT linear. Integrate this factor
at multiple Gauss-Legendre points instead of evaluating only at the
midpoint. Analytically differentiate the resulting discrete ACTION to
obtain momenta, mass matrix, position forces and convective derivatives.

This is a controlled discretization/aliasing test, *not* a fully
continuum-validated Nambu--Goto solver. Equal-arclength gauge and its
limitations are inherited from EquidistributedNGString.
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

from .equidistributed_ng import EquidistributedNGString
from .dynamics import C_CM_PER_PS,HBARC_GEV_CM
from .field_dynamics import MAGNETIC_FORCE_GEV_PER_PS_T

I3=np.eye(3)


class GaussQuadratureNGString(EquidistributedNGString):
    def __init__(self,params,quadrature_order=2):
        super().__init__(params)
        if quadrature_order not in (1,2,3):
            raise ValueError('quadrature_order must be 1, 2 or 3')
        self.quadrature_order=quadrature_order
        q,w=leggauss(quadrature_order)
        self.gauss_u=(q+1.)/2.
        self.gauss_w=w/2.

    def action_quantities(self,state,derivatives=False):
        x,v=self._unpack(state)
        n=self.n
        v2=np.einsum('ij,ij->i',v,v)
        if np.max(v2)>=1:
            raise ValueError('superluminal discrete string vertex')
        gamma=1/np.sqrt(1-v2)
        canonical=self.masses[:,None]*gamma[:,None]*v
        energy=float(np.sum(self.masses*gamma))
        matrix=np.zeros((3*n,3*n)) if derivatives else None
        force=np.zeros((n,3)) if derivatives else None
        convection=np.zeros((n,3)) if derivatives else None
        if derivatives:
            for i in range(n):
                sl=slice(3*i,3*i+3)
                matrix[sl,sl]+=self.masses[i]*gamma[i]*(
                    I3+gamma[i]**2*np.outer(v[i],v[i]))
        for i in range(n-1):
            d=x[i+1]-x[i]
            ell=math.sqrt(float(d@d)+self.params.regularizer_cm**2)
            tangent=d/ell
            d_dot=C_CM_PER_PS*(v[i+1]-v[i])
            ell_dot=float(tangent@d_dot)
            tan_dot=(d_dot-ell_dot*tangent)/ell
            left=slice(3*i,3*i+3)
            right=slice(3*(i+1),3*(i+1)+3)
            for u,weight in zip(self.gauss_u,self.gauss_w):
                f0=1.-u
                f1=u
                vel=f0*v[i]+f1*v[i+1]
                par=float(vel@tangent)
                perp=vel-par*tangent
                h2=1-float(perp@perp)
                if h2<=0:
                    raise ValueError('non-timelike spatial quadrature point')
                h=math.sqrt(h2)
                coeff=self.tension*weight
                ps=coeff*ell*perp/h
                canonical[i]+=f0*ps
                canonical[i+1]+=f1*ps
                energy+=coeff*ell/h
                if not derivatives:
                    continue
                Q=I3-np.outer(tangent,tangent)
                block=coeff*ell*(Q/h+np.outer(perp,perp)/h**3)
                matrix[left,left]+=f0*f0*block
                matrix[left,right]+=f0*f1*block
                matrix[right,left]+=f0*f1*block
                matrix[right,right]+=f1*f1*block
                # Variation of L w.r.t. endpoints x_i and x_(i+1).
                ff=C_CM_PER_PS*coeff*(h*tangent+(par/h)*perp)
                force[i]+=ff
                force[i+1]-=ff
                par_dot=float(vel@tan_dot)
                perp_dot=-par_dot*tangent-par*tan_dot
                h_dot=par*par_dot/h
                dot=coeff*((ell_dot/h-ell*h_dot/h**2)*perp+
                           ell/h*perp_dot)
                convection[i]+=f0*dot
                convection[i+1]+=f1*dot
        result={
            'energy_gev':float(energy),
            'momentum_gev':np.sum(canonical,axis=0),
            'angular_momentum_over_hbar':
                np.sum(np.cross(x,canonical),axis=0)/HBARC_GEV_CM,
            'canonical_vertex_momenta_gev':canonical
        }
        if derivatives:
            B=np.array([0.,0.,self.params.b_tesla])
            for endpoint in (0,n-1):
                force[endpoint]+=self.charge[endpoint]*(
                    MAGNETIC_FORCE_GEV_PER_PS_T*np.cross(v[endpoint],B))
            result.update({
                'mass_matrix':matrix,
                'forces_gev_per_ps':force,
                'convective_gev_per_ps':convection
            })
        return result
