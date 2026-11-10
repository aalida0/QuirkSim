"""Fail-closed element-quality events for experimental finite-element worldsheets.

Quality thresholds are *engineering* diagnostics, not physical distances.
They stop a solve at the FIRST underresolved element. Low drift in global
Noether charges does not license propagating through a degenerate mesh.
"""
import numpy as np


class MeshQualityViolation(RuntimeError):
    pass


def mesh_metrics(state,nodes,regulator_cm):
    pos=np.asarray(state[:3*nodes],float).reshape(nodes,3)
    lengths=np.linalg.norm(np.diff(pos,axis=0),axis=1)
    min_length=float(lengths.min())
    max_length=float(lengths.max())
    return {
        'min_element_cm':min_length,'max_element_cm':max_length,
        'max_to_min_ratio':max_length/max(min_length,regulator_cm)
    }


def is_mesh_acceptable(state,nodes,regulator_cm,
                       minimum_factor=100.,maximum_ratio=50.):
    m=mesh_metrics(state,nodes,regulator_cm)
    return bool(m['min_element_cm']>=minimum_factor*regulator_cm
                and m['max_to_min_ratio']<=maximum_ratio)


def make_mesh_stop_events(nodes,regulator_cm,
                          minimum_factor=100.,maximum_ratio=50.):
    """Pass both events to solve_ivp in either velocity or rapidity variables."""
    def short_element(_t,y):
        m=mesh_metrics(y,nodes,regulator_cm)
        return m['min_element_cm']-minimum_factor*regulator_cm

    def imbalance(_t,y):
        m=mesh_metrics(y,nodes,regulator_cm)
        return maximum_ratio-m['max_to_min_ratio']

    for fn in (short_element,imbalance):
        fn.terminal=True
        fn.direction=-1
    return [short_element,imbalance]


def require_good_mesh(state,nodes,regulator_cm):
    if not is_mesh_acceptable(state,nodes,regulator_cm):
        raise MeshQualityViolation(
            f'initial mesh not acceptable: {mesh_metrics(state,nodes,regulator_cm)}')


def assert_no_mesh_stop(solution):
    kinds=('short_element','element_imbalance')
    for name, hits in zip(kinds,solution.t_events[-2:]):
        if len(hits):
            raise MeshQualityViolation(
                f'worldsheet reached {name} gate at t={float(hits[0]):.9g} ps; '
                'no physical results are valid beyond this time')


def total_string_arclength_cm(state,nodes):
    """Sum segment lengths irrespective of interior-label spacing.

    This quantity is the physical spatial arclength of the represented
    piecewise-linear string at a single lab time, NOT an annihilation radius.
    """
    pos=np.asarray(state[:3*nodes],dtype=float).reshape(nodes,3)
    return float(np.sum(np.linalg.norm(np.diff(pos,axis=0),axis=1)))


def make_total_arclength_stop_event(nodes,cutoff_cm):
    """Stop contracting strings at a common TOTAL arclength.

    Unlike a per-segment threshold this is directly comparable for
    different numbers of vertices, even with imperfect gauge spacing.
    It does not define a physical continuation through a pinch.
    """
    if not np.isfinite(cutoff_cm) or cutoff_cm<=0:
        raise ValueError('total arclength cutoff must be positive and finite')

    def total_arclength_pinch(_t,state):
        return total_string_arclength_cm(state,nodes)-cutoff_cm

    total_arclength_pinch.terminal=True
    total_arclength_pinch.direction=-1
    return total_arclength_pinch
