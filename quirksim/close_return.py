"""Conservative post-run close-return censor for asynchronous endpoint samples.

In exact arithmetic, the bound assumes continuous worldlines with speed <= c.
It does not establish a crossing, a quantum contact, or annihilation.
"""
import numpy as np

C_MM_PER_PS = 0.299792458


def speed_cone_window(t0, x0, t1, x1, max_grid_spacing_ps=3e-5,
                      roundoff_margin_mm=1e-9):
    """Bound the minimum separation over overlapping recorded histories.

    At t, each observed position confines its endpoint within c*|t-t_obs|.
    Maximize four triangle-inequality lower bounds, using the adjacent sample
    on each side for each endpoint. The actual distance is 2c-Lipschitz; a
    uniform grid with spacing h therefore loses at most c*h between queries.
    This assumes continuous subluminal histories, with a heuristic roundoff
    guard. It neither interpolates a trajectory nor certifies transport error.
    """
    if (not np.isfinite(max_grid_spacing_ps) or max_grid_spacing_ps <= 0
            or not np.isfinite(roundoff_margin_mm) or roundoff_margin_mm < 0):
        raise ValueError('positive spacing and nonnegative roundoff guard required')
    a, xa = _endpoint(t0,x0); b, xb = _endpoint(t1,x1)
    lo=max(a[0],b[0]); hi=min(a[-1],b[-1])
    if lo >= hi: raise ValueError('no overlapping time span')
    cells=int(np.ceil((hi-lo)/max_grid_spacing_ps))
    if cells > 10000000: raise ValueError('window too large: split into overlapping subwindows')
    spacing=(hi-lo)/cells
    grid=np.unique(np.r_[np.linspace(lo,hi,cells+1), a[(a>=lo)&(a<=hi)], b[(b>=lo)&(b<=hi)]])
    minimum_lower=np.inf; minimum_upper=np.inf; witness=None
    for offset in range(0,len(grid),100000):
        g=grid[offset:offset+100000]
        ia=np.clip(np.searchsorted(a,g),1,len(a)-1)
        ib=np.clip(np.searchsorted(b,g),1,len(b)-1)
        lower=np.zeros(len(g)); upper=np.full(len(g),np.inf)
        for da in (-1,0):
            for db in (-1,0):
                j,k=ia+da,ib+db
                distance=np.linalg.norm(xa[j]-xb[k],axis=1)
                travel=C_MM_PER_PS*(np.abs(g-a[j])+np.abs(g-b[k]))
                lower=np.maximum(lower,distance-travel)
                upper=np.minimum(upper,distance+travel)
        minimum_lower=min(minimum_lower,float(lower.min()))
        k=int(np.argmin(upper))
        if upper[k] < minimum_upper:
            minimum_upper=float(upper[k]); witness=float(g[k])
    return {'common_window_ps':[float(lo),float(hi)],
            'minimum_separation_lower_bound_mm':max(0.,minimum_lower-C_MM_PER_PS*spacing-roundoff_margin_mm),
            'minimum_separation_upper_bound_mm':minimum_upper+roundoff_margin_mm,
            'upper_bound_witness_time_ps':witness,'grid_spacing_ps':spacing,
            'roundoff_margin_mm':roundoff_margin_mm,
            'meaning':'SUBLUMINAL_RECORDED_WORLDLINES_ONLY_NOT_QUANTUM_SURVIVAL',
            'arithmetic':'FLOATING_POINT_MARGIN_NOT_FORMALLY_CERTIFIED'}


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
