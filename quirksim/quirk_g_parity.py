"""Twisted charge conjugation for the minimal (2_IC,3_EW,Y=0) Dirac model.

G acts on the quirk through epsilon_IC tensor S_EW times charge conjugation;
all gauge fields can be held fixed for the two-color theory. No new breaking
operators, hypercharge or chiral heavy-current couplings are included.
"""
import numpy as np
from .production_spin import PAULI


def generators():
    plus=np.sqrt(2)*np.array([[0,1,0],[0,0,1],[0,0,0]],complex)
    return np.array([(plus+plus.T)/2,(plus-plus.T)/(2j),np.diag([1.,0.,-1.])])


def invariant_checks():
    ew=generators();s=np.array([[0,0,1],[0,-1,0],[1,0,0]],complex)
    epsilon=np.array([[0,1],[-1,0]],complex);a=np.kron(epsilon,s)
    gauge=[np.kron(t/2,np.eye(3)) for t in PAULI]+[np.kron(np.eye(2),t) for t in ew]
    return {'gauge_conjugation_residual':float(max(np.max(abs(a@t.conj()@a.conj().T+t)) for t in gauge)),
        'double_conjugation_residual':float(np.max(abs(a@a.conj()+np.eye(6)))),
        'mass_splitting_residual':float(np.max(abs(s@np.diag([1.,0.,1.])@s.conj().T-np.diag([1.,0.,1.])))),
        'scope':'minimal gauge kinetic and symmetric split-mass algebra; excludes G-breaking higher-dimension operators'}


def neutral_pair_g_parity(orbital_l,spin):
    if type(orbital_l) is not int or orbital_l<0 or spin not in (0,1):raise ValueError('integer L>=0 and S=0/1 required')
    # Flavor reversal leaves the neutral particle-antiparticle channel fixed.
    return (-1)**(orbital_l+spin)


def neutral_capture_channels(source_g=1):
    if source_g not in (-1,1):raise ValueError('G must be +/-1')
    return {name:{'L':ell,'S':spin,'G':neutral_pair_g_parity(ell,spin),
            'allowed_from_source_with_G_even_radiation':neutral_pair_g_parity(ell,spin)==source_g}
            for name,ell,spin in [('1S0',0,0),('3S1',0,1),('3P0',1,1),('3P1',1,1),('3P2',1,1)]}
