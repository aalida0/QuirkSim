"""Exact 1+1D radial Nambu--Goto straight-string massive-endpoint benchmark.

Sector: two identical massive endpoints, equal opposite momenta, zero
transverse degrees, B=0, no weak decay or radiation. In the strictly radial
target-time ansatz, the string Lagrangian is -T * |x_R-x_L|, with
T = Lambda_IC**2/(hbar*c). The endpoints satisfy

  d p_R / dt = -sign(x_R-x_L) * Lambda_IC**2/hbar,
  d p_L / dt = -d p_R / dt.

Within each interval between zero-length crossings the force is CONSTANT.
Relativistic endpoint positions and proper times have exact primitives.
We continue through crossings as a signed 1D ideal-string model; this is
an explicit *crossing convention*, NOT a unique 3D reconnection/annihilation
prescription nor a full curved Nambu--Goto worldsheet evolution.

Use this as an independent multi-oscillation reference for convergence
and gauge pathology diagnostics. Do not use as recoil or LHC prediction.
"""
from dataclasses import dataclass
import math
import numpy as np

from .physics import benchmark, HBAR_GEV_PS
from .dynamics import C_CM_PER_PS, HBARC_GEV_CM


@dataclass(frozen=True)
class RadialNGParameters:
    mass_neutral_gev: float = 225.0
    lambda_ev: float = 774.0
    initial_beta: float = 0.4
    initial_separation_cm: float = 0.1

    def __post_init__(self):
        if not (self.mass_neutral_gev>0 and self.lambda_ev>0
                and 0<=self.initial_beta<1
                and self.initial_separation_cm>=0):
            raise ValueError('invalid radial reference benchmark')


class ExactRadialNG:
    """Piecewise analytic radial two-endpoint/string Hamiltonian reference."""

    def __init__(self,p=None):
        if p is None:
            p=RadialNGParameters()
        self.params=p
        self.mass=benchmark('triplet',p.mass_neutral_gev)['charged_mass_gev']
        self.tension=(p.lambda_ev*1e-9)**2/HBARC_GEV_CM
        self.force=(p.lambda_ev*1e-9)**2/HBAR_GEV_PS
        self.p0=self.mass*p.initial_beta/math.sqrt(1-p.initial_beta**2)
        self.energy_total=2*math.hypot(self.mass,self.p0)+(
            self.tension*p.initial_separation_cm)
        e_cross=self.energy_total/2
        # Factorized for small initial string energies: avoids cancellation.
        self.p_cross=math.sqrt(max(0,(e_cross-self.mass)*(e_cross+self.mass)))
        if self.p_cross==0:
            raise ValueError('degenerate completely static/zero-string benchmark')
        self.first_crossing_ps=(self.p0+self.p_cross)/self.force
        self.crossing_interval_ps=2*self.p_cross/self.force
        self.max_separation_cm=(self.energy_total-2*self.mass)/self.tension
        self.first_segment_proper_ps=(self.mass/self.force)*(
            math.asinh(self.p0/self.mass)+math.asinh(self.p_cross/self.mass))
        self.subsequent_segment_proper_ps=2*self.mass/self.force*(
            math.asinh(self.p_cross/self.mass))

    def _branch(self,t):
        if not math.isfinite(t) or t<0:
            raise ValueError('time must be finite and nonnegative')
        if t < self.first_crossing_ps:
            crossings=0
            phase=t
            momentum=self.p0-self.force*t
            sign=1
        else:
            elapsed=t-self.first_crossing_ps
            crossings=1+math.floor(elapsed/self.crossing_interval_ps)
            phase=elapsed-(crossings-1)*self.crossing_interval_ps
            # Crossing 1: right endpoint moving left; force flips outward.
            # Crossing 2: right endpoint moving right; force flips inward.
            if crossings%2:
                momentum=-self.p_cross+self.force*phase
                sign=-1
            else:
                momentum=self.p_cross-self.force*phase
                sign=1
        return crossings,phase,momentum,sign

    def state(self,t):
        """Return exact endpoint positions/momenta, string energy, proper time.

        The unsigned separation is always nonnegative. 'crossings' counts
        endpoint worldline intersections since t=0; no annihilation is assumed.
        """
        crossing,phase,momentum,sign=self._branch(t)
        energy=math.hypot(self.mass,momentum)
        length=max(0,(self.energy_total-2*energy)/self.tension)
        signed_length=sign*length
        if crossing==0:
            tau=(self.mass/self.force)*(
                math.asinh(self.p0/self.mass)-math.asinh(momentum/self.mass))
        else:
            if crossing%2:
                tau_partial=(self.mass/self.force)*(
                    math.asinh(momentum/self.mass)+
                    math.asinh(self.p_cross/self.mass))
            else:
                tau_partial=(self.mass/self.force)*(
                    math.asinh(self.p_cross/self.mass)-
                    math.asinh(momentum/self.mass))
            tau=self.first_segment_proper_ps+(
                (crossing-1)*self.subsequent_segment_proper_ps+tau_partial)
        result={
            'time_ps':float(t),
            'right_position_cm':signed_length/2,
            'left_position_cm':-signed_length/2,
            'right_momentum_gev':momentum,
            'left_momentum_gev':-momentum,
            'right_beta':momentum/energy,
            'left_beta':-momentum/energy,
            'proper_time_each_ps':tau,
            'endpoint_energy_each_gev':energy,
            'string_energy_gev':self.tension*length,
            'total_energy_gev':2*energy+self.tension*length,
            'total_linear_momentum_gev':0.,
            'total_angular_momentum_hbar':0.,
            'crossings':crossing,
            'unsigned_separation_cm':length,
            'signed_relative_position_cm':signed_length,
            'status':'EXACT_REDUCED_RADIAL_NG_WITH_SIGNED_CROSSING_CONVENTION'
        }
        return result

    def gauge_slice(self,t,fractions):
        """Exact radial worldsheet snapshot with an arbitrary monotone labelling.

        In the strictly 1D sector, longitudinal string motion is a gauge
        direction. Any fixed monotone u in [0,1] parameterizes the *same*
        straight radial string at the given physical time. Each material
        worldsheet vertex has velocity xdot(u)=(1-u)xdot_L+u xdot_R,
        with NO artificial vertex rest masses. Consequently the string's
        physical canonical longitudinal momentum vanishes identically,
        regardless of label spacing, and its total energy is T*|r|.

        This does NOT prescribe remapping or conserved canonical momenta
        of a curved, transversely excited 3D worldsheet.
        """
        u=np.asarray(fractions,dtype=float)
        if (u.ndim!=1 or len(u)<2 or not np.all(np.isfinite(u))
                or abs(float(u[0]))>1e-13 or abs(float(u[-1])-1)>1e-13
                or np.any(np.diff(u)<=0)):
            raise ValueError('strictly increasing gauge fractions from 0 to 1 needed')
        state=self.state(t)
        right=state['right_position_cm']
        left=state['left_position_cm']
        br=state['right_beta']
        bl=state['left_beta']
        positions=(1-u)*left+u*right
        velocities=(1-u)*bl+u*br
        lengths=np.abs(np.diff(positions))
        string_energy=float(self.tension*np.sum(lengths))
        end_momentum=(state['left_momentum_gev']+
                      state['right_momentum_gev'])
        total_energy=2*state['endpoint_energy_each_gev']+string_energy
        return {
            'time_ps':t,
            'fractions':u.tolist(),
            'vertex_positions_cm':positions.tolist(),
            'vertex_betas':velocities.tolist(),
            'endpoint_momentum_gev':end_momentum,
            'string_longitudinal_momentum_gev':0.0,
            'string_energy_gev':string_energy,
            'total_energy_gev':total_energy,
            'max_element_ratio':float(np.max(lengths)/np.min(lengths))
                if np.min(lengths)>0 else None,
            'gauge_status':'EXACT_RADIAL_REPARAMETRIZATION_ONLY',
            'limits':'not a conservative remesh of a moving 3D worldsheet'
        }

    def summary(self,duration_ps,samples=201):
        if samples<2:
            raise ValueError('samples must be >= 2')
        times=np.linspace(0,duration_ps,samples)
        values=[self.state(float(t)) for t in times]
        deviations=[abs(z['total_energy_gev']-self.energy_total) for z in values]
        return {
            'reference':'exact 1D signed radial massive-endpoint NG sector',
            'mass_gev':self.mass,
            'lambda_ev':self.params.lambda_ev,
            'initial_separation_cm':self.params.initial_separation_cm,
            'first_crossing_ps':self.first_crossing_ps,
            'crossing_interval_ps':self.crossing_interval_ps,
            'duration_ps':duration_ps,
            'crossings_at_end':values[-1]['crossings'],
            'max_relative_energy_residual':max(deviations)/self.energy_total,
            'maximum_separation_cm':self.max_separation_cm,
            'proper_time_each_endpoint_at_end_ps':
                values[-1]['proper_time_each_ps'],
            'final_unsigned_separation_cm':
                values[-1]['unsigned_separation_cm'],
            'limits':[
                'radial 1D symmetry; no curvature, transverse modes or radiation',
                'signed crossing continuation is a convention not annihilation',
                'not a moving 3D worldsheet benchmark or physical pp initial state']
        }
