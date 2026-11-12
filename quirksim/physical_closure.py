"""Frozen neutral-vacuum phenomenological closure, not a unique quantum rate."""
import math
from .physics import HBAR_GEV_PS


def compact_subspace_interval(capture_probability,capture_age_ps,widths_gev):
    """Any density supported on a specified diagonal compact decay subspace.

    Requires the supplied widths to cover the support; it does not determine
    capture occupation or include unlisted/metastable states.
    """
    widths=tuple(widths_gev)
    if not widths or any(not math.isfinite(w) or w<=0 for w in widths):
        raise ValueError('positive finite width for every admitted compact sector required')
    low=capture_annihilation_interval(capture_probability,capture_age_ps,min(widths),1.)['lower']
    high=capture_annihilation_interval(capture_probability,capture_age_ps,max(widths),1.)['upper']
    return {'lower':low,'upper':high,'scope':'SPECIFIED_COMPACT_SUBSPACE_OCCUPATION_ENVELOPE'}


def capture_annihilation_interval(capture_probability,capture_age_ps,width_gev,
                                  singlet_fraction=None):
    """Keep uncalculated spin-sector annihilation unknown, rather than zero.

    Capture is assumed to populate the specified compact-state family. No
    pre-capture annihilation is included in this conditional interval.
    """
    if (not all(math.isfinite(v) for v in (capture_probability,capture_age_ps,width_gev))
        or not 0<=capture_probability<=1 or capture_age_ps<0 or width_gev<=0
        or singlet_fraction is not None and (not math.isfinite(singlet_fraction) or not 0<=singlet_fraction<=1)):
        raise ValueError('invalid capture/spin/width inputs')
    decay=-math.expm1(-capture_age_ps*width_gev/HBAR_GEV_PS)
    if capture_age_ps==0:
        return {'lower':0.,'upper':0.,'spin_status':'ZERO_POST_CAPTURE_AGE'}
    if singlet_fraction is None:
        return {'lower':0.,'upper':capture_probability,'spin_status':'UNSPECIFIED'}
    known=capture_probability*singlet_fraction*decay
    return {'lower':known,'upper':known+capture_probability*(1-singlet_fraction),
            'spin_status':'SPECIFIED_SINGLET_FRACTION_OTHER_SECTOR_UNCALCULATED'}


def centrifugal_action_floor(mass_gev,excess_energy_gev,angular,
                              attractive_coulomb_ceiling=.2,string_tension_gev2=(774e-9)**2):
    """Action lower bound for a specified large-L semiclassical Hamiltonian.

    H=2 sqrt(m^2+p_r^2+L^2/r^2)+sigma*r-kappa/r, 0<=kappa<=ceiling.
    On 1/m <= r <= L/H and L>=2*kappa, |p_r|>=sqrt(7)*L/(4*r).
    The exponential is a WKB diagnostic; it is not a rigorous probability
    bound for an unspecified wavepacket or a measured partial-wave spectrum.
    """
    if not all(math.isfinite(x) for x in (mass_gev,excess_energy_gev,angular,attractive_coulomb_ceiling,string_tension_gev2)) or min(mass_gev,excess_energy_gev)<=0 or min(angular,attractive_coulomb_ceiling,string_tension_gev2)<0:
        raise ValueError('invalid Hamiltonian inputs')
    ratio=mass_gev*angular/(2*mass_gev+excess_energy_gev)
    if angular<2*attractive_coulomb_ceiling or ratio<=1 or string_tension_gev2*angular>=(2*mass_gev+excess_energy_gev)**2:
        return {'status':'NOT_IN_LARGE_L_DOMAIN','action_lower_bound':None,'log10_WKB_exponential_ceiling':None}
    action=math.sqrt(7)*angular*math.log(ratio)/4
    return {'status':'SPECIFIED_CENTRAL_WKB_DIAGNOSTIC','action_lower_bound':action,
            'log10_WKB_exponential_ceiling':-2*action/math.log(10),
            'attractive_coulomb_ceiling':attractive_coulomb_ceiling}
