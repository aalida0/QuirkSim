"""Four-momentum conserving weak-decay kinematics of unpolarized Dirac quirks.

Charged vectorlike quirk -> neutral quirk + pion or charged lepton + neutrino.
Leading four-Fermi vector heavy current times V-A SM lepton current.

MOMENTUM SHAPES: Monte Carlo rejection sampling of the exact tree-level
three-body decay matrix-element shape for an unpolarized vectorlike current.
TOTAL RATES: supplied by the separate LO width formulas in physics.py.
No finite-W, EW NLO, polarization, bound-state or IC corrections.
"""
import math
import numpy as np
from .physics import MPI, ME, MMU


def boost(e, momentum, beta):
    """Lorentz boost of four-vector (E, p) by velocity beta, natural units."""
    momentum = np.asarray(momentum, dtype=float)
    beta = np.asarray(beta, dtype=float)
    b2 = float(beta @ beta)
    if b2 == 0:
        return np.array([e, *momentum], dtype=float)
    if b2 >= 1:
        raise ValueError('superluminal Lorentz boost')
    gamma = 1 / math.sqrt(1 - b2)
    bp = float(beta @ momentum)
    return np.array([gamma * (e + bp), *(momentum + ((gamma - 1) * bp / b2 + gamma * e) * beta])


def direction(rng):
    costheta = rng.uniform(-1, 1)
    phi = rng.uniform(0, 2 * math.pi)
    s = math.sqrt(1 - costheta**2)
    return np.array([s*math.cos(phi), s*math.sin(phi), costheta])


def scalar(a, b):
    return float(a[0]*b[0] - a[1:] @ b[1:])


def vector_current_weight(parent, daughter, lepton, neutrino):
    """Unpolarized vector x V-A spin-summed shape, overall constant removed.

    Symmetric tensor contraction proportional to:
    (p0.pL)(pP.pNu) + (p0.pNu)(pP.pL) - M*m*(pL.pNu).
    Antisymmetric terms drop for an unpolarized Dirac vector current.
    """
    m = math.sqrt(max(scalar(daughter, daughter), 0))
    M = math.sqrt(max(scalar(parent, parent), 0))
    return (scalar(daughter, lepton) * scalar(parent, neutrino)
            + scalar(daughter, neutrino) * scalar(parent, lepton)
            - M*m*scalar(lepton, neutrino))


def pion_decay_rest(M, m, rng):
    if M <= m + MPI:
        raise ValueError('pion channel not kinematically accessible')
    k = math.sqrt((M**2-(m+MPI)**2)*(M**2-(m-MPI)**2))/(2*M)
    n = direction(rng)
    return {
        'daughter': np.array([math.hypot(m,k), *(-k*n)]),
        'pion': np.array([math.hypot(MPI,k), *(k*n)]),
    }


def leptonic_decay_rest(M, m, ml, rng, max_tries=20000):
    """Sample exact three-body unpolarized V x V-A shape in parent rest frame.

    Generate q^2 uniformly, q direction and lepton angle uniformly, then
    accept with Jacobian p(q) * (1-ml^2/q^2) times the tensor contraction.
    An analytic conservative upper bound prevents adaptive-rejection bias.
    """
    if not 0 <= ml < M-m:
        raise ValueError('leptonic channel not accessible')
    delta = M - m
    parent = np.array([M, 0., 0., 0.])
    bound = 2 * M * (M + delta) * delta**3
    lower, upper = ml**2, delta**2
    for _ in range(max_tries):
        s = rng.uniform(lower, upper)
        sqrt_s = math.sqrt(s)
        if sqrt_s == 0:
            continue
        root_lambda = math.sqrt(max(0., (delta**2-s)*((M+m)**2-s)))
        pq = root_lambda / (2*M)
        eq = (M*M-m*m+s)/(2*M)
        ed = (M*M+m*m-s)/(2*M)
        n = direction(rng)
        k = (s-ml*ml)/(2*sqrt_s)
        el = (s+ml*ml)/(2*sqrt_s)
        lstar = k * direction(rng)
        lb = boost(el, lstar, pq*n/eq)
        nb = boost(k, -lstar, pq*n/eq)
        d = np.array([ed, *(-pq*n)])
        matrix = vector_current_weight(parent, d, lb, nb)
        if matrix < -1e-7:
            raise RuntimeError(f'negative matrix-element weight: {matrix}')
        weight = pq * (1-ml*ml/s) * max(0.,matrix)
        if weight > bound*(1+1e-10):
            raise RuntimeError('rejection envelope violated')
        if rng.random() * bound < weight:
            return {'daughter': d, 'charged_lepton': lb, 'neutrino': nb}
    raise RuntimeError('leptonic phase-space rejection failed; check mass gap')


def sample_decay(channel, parent_momentum, parent_mass, daughter_mass, rng):
    """Return daughter 4-vector and visible/invisible products in lab frame."""
    parent_momentum = np.asarray(parent_momentum, float)
    parent_energy = math.hypot(parent_mass, np.linalg.norm(parent_momentum))
    beta = parent_momentum/parent_energy
    if channel == 'pion':
        rest = pion_decay_rest(parent_mass, daughter_mass, rng)
    elif channel == 'electron_neutrino':
        rest = leptonic_decay_rest(parent_mass, daughter_mass, ME, rng)
    elif channel == 'muon_neutrino':
        rest = leptonic_decay_rest(parent_mass, daughter_mass, MMU, rng)
    else:
        raise ValueError(f'unknown decay channel: {channel}')
    lab = {name: boost(v[0], v[1:], beta) for name,v in rest.items()}
    initial = np.array([parent_energy, *parent_momentum])
    final = sum(lab.values(), np.zeros(4))
    if not np.allclose(initial, final, rtol=0, atol=3e-9):
        raise RuntimeError(f'weak-decay four-momentum violation: {initial-final}')
    return lab
