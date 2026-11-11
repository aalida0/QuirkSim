"""Conservative post-run close-return censor for asynchronous endpoint samples.

In exact arithmetic, the bound assumes continuous worldlines with speed <= c.
It does not establish a crossing, a quantum contact, or annihilation.
"""
import numpy as np

C_MM_PER_PS = 0.299792458


def _bound(a, ta, b, tb, p, tp):
    dt = tb - ta
    travel = C_MM_PER_PS * np.maximum(np.abs(ta-tp), np.abs(tb-tp))
    return np.maximum(0., np.maximum(np.linalg.norm(a-p, axis=1),
                                     np.linalg.norm(b-p, axis=1))
                      - C_MM_PER_PS*dt - travel)


def _endpoint(times, positions):
    t, x = np.asarray(times, float), np.asarray(positions, float)
    if t.ndim != 1 or x.shape != (len(t), 3) or len(t) < 2:
        raise ValueError('each endpoint needs at least two 3D samples')
    if not np.isfinite(t).all() or not np.isfinite(x).all() or np.any(np.diff(t) < 0):
        raise ValueError('nonfinite or time-reversed endpoint samples')
    # Handoffs can produce two samples at the same time. They must coincide.
    dt = np.diff(t)
    dx = np.linalg.norm(np.diff(x, axis=0), axis=1)
    if np.any(dx > C_MM_PER_PS*dt + 1e-10):
        raise ValueError('discontinuous or superluminal sampled endpoint')
    keep = np.r_[dt > 0, True]
    t, x = t[keep], x[keep]
    if len(t) < 2: raise ValueError('endpoint samples have no time span')
    return t, x


def monitor(t0, x0, t1, x1, *, start_ps, arm_radius_mm, match_radius_mm,
            roundoff_margin_mm=1e-9):
    """Find earliest unproven interval after the pair is demonstrably separated.

    Both lists must cover the requested start and represent the entire endpoint
    history, including weak-decay handoffs. Output concerns recorded transport
    only. The margin is heuristic, not a formal floating-point error bound.
    """
    if not all(np.isfinite(v) for v in (start_ps, arm_radius_mm, match_radius_mm, roundoff_margin_mm)):
        raise ValueError('nonfinite monitor configuration')
    if not (arm_radius_mm > match_radius_mm > 0 and roundoff_margin_mm >= 0):
        raise ValueError('require arm radius > matching radius > 0 and nonnegative margin')
    a, xa = _endpoint(t0, x0); b, xb = _endpoint(t1, x1)
    begin=max(a[0], b[0], start_ps); end=min(a[-1], b[-1])
    result={'status':'INSUFFICIENT_OVERLAPPING_HISTORY', 'requested_start_ps':float(start_ps),
            'arm_radius_mm':float(arm_radius_mm), 'match_radius_mm':float(match_radius_mm),
            'roundoff_margin_mm':float(roundoff_margin_mm),
            'meaning':'NUMERICAL_TRANSPORT_CENSOR_NOT_CONTACT_OR_ANNIHILATION',
            'arithmetic':'FLOATING_POINT_MARGIN_NOT_FORMALLY_CERTIFIED'}
    if begin >= end: return result
    # Common time cells: each is enclosed by one step of each endpoint.
    boundaries=np.unique(np.r_[begin, a[(a>begin)&(a<end)], b[(b>begin)&(b<end)], end])
    left, right=boundaries[:-1], boundaries[1:]
    ia=np.searchsorted(a,left,side='right')-1; ib=np.searchsorted(b,left,side='right')-1
    # Use both partner observations and both endpoint orderings, keeping the
    # strongest lower bound. Each evaluated bound covers a whole source step.
    bounds=[]
    for j in (0,1):
        bounds.append(_bound(xa[ia],a[ia],xa[ia+1],a[ia+1],xb[ib+j],b[ib+j]))
        bounds.append(_bound(xb[ib],b[ib],xb[ib+1],b[ib+1],xa[ia+j],a[ia+j]))
    lower=np.maximum.reduce(bounds)-roundoff_margin_mm
    arm=np.flatnonzero(lower>arm_radius_mm)
    result.update(common_start_ps=float(begin), common_end_ps=float(end), cells=int(len(left)))
    if not len(arm):
        result['status']='NEVER_ARMED'; return result
    first=int(arm[0]); result['armed_at_ps']=float(left[first])
    uncertain=np.flatnonzero(lower[first:]<=match_radius_mm)
    if len(uncertain):
        k=first+int(uncertain[0]);result.update(status='CENSORED_UNRESOLVED_INTERVAL',
            censor_interval_ps=[float(left[k]),float(right[k])],lower_bound_mm=float(lower[k]))
    else:
        result.update(status='OUTSIDE_RADIUS_UNTIL_COMMON_END',
                      lower_bound_mm=float(np.min(lower[first:])))
    return result
