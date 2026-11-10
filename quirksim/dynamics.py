"""Relativistic two-endpoint straight-string *proxy* with weak pion emission.

Units: GeV for momenta/masses, cm for positions, ps for laboratory time.
NOT a Nambu--Goto worldsheet, radiation or annihilation simulation.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.integrate import solve_ivp
from .physics import benchmark, MPI, HBAR_GEV_PS

C_CM_PER_PS = 0.0299792458


@dataclass(frozen=True)
class Parameters:
    representation: str = "triplet"
    mass_neutral_gev: float = 225.0
    lambda_ev: float = 774.0
    beta: float = 0.4
    impact_cm: float = 0.0
    regulator_cm: float = 1e-7
    duration_ps: float = 500.0
    seed: int = 1
    mode: str = "recoil"

    def __post_init__(self):
        if self.representation not in ("triplet", "doublet"):
            raise ValueError("unknown representation")
        if self.mode not in ("stable", "no_recoil", "recoil"):
            raise ValueError("unknown mode")
        if self.mass_neutral_gev <= 0 or self.lambda_ev <= 0 or self.duration_ps <= 0:
            raise ValueError("nonpositive physical parameter")
        if not 0 < self.beta < 1 or self.impact_cm < 0 or self.regulator_cm <= 0:
            raise ValueError("invalid initial conditions")


def boost(e, p, beta):
    beta = np.asarray(beta, float)
    p = np.asarray(p, float)
    b2 = float(beta @ beta)
    if b2 == 0:
        return float(e), p.copy()
    if b2 >= 1:
        raise ValueError("invalid beta")
    gamma = 1.0 / math.sqrt(1.0-b2)
    bp = float(beta @ p)
    return gamma*(e+bp), p + ((gamma-1)*bp/b2 + gamma*e)*beta


def pion_two_body(momentum, parent_mass, daughter_mass, rng):
    """Unpolarized isotropic two-body decay in parent rest frame, boosted to lab."""
    momentum = np.asarray(momentum, float)
    parent_e = math.hypot(parent_mass, np.linalg.norm(momentum))
    kallen = ((parent_mass**2 - (daughter_mass+MPI)**2) *
              (parent_mass**2 - (daughter_mass-MPI)**2))
    if kallen <= 0:
        raise ValueError("pion decay is kinematically closed")
    pstar = math.sqrt(kallen)/(2*parent_mass)
    costheta = rng.uniform(-1,1)
    phi = rng.uniform(0,2*math.pi)
    n = np.array([math.sqrt(1-costheta**2)*math.cos(phi),
                  math.sqrt(1-costheta**2)*math.sin(phi), costheta])
    beta = momentum/parent_e
    daughter = boost(math.hypot(daughter_mass,pstar),-pstar*n,beta)
    pion = boost(math.hypot(MPI,pstar),pstar*n,beta)
    if not np.allclose(daughter[1]+pion[1],momentum,atol=1e-9):
        raise RuntimeError("decay violates momentum conservation")
    if abs(daughter[0]+pion[0]-parent_e)>1e-9:
        raise RuntimeError("decay violates energy conservation")
    return daughter[1], pion


def total_energy(state, masses, sigma, regulator_cm):
    r = np.linalg.norm(state[:3]-state[3:6])
    string = sigma/(C_CM_PER_PS*HBAR_GEV_PS)*(math.hypot(r,regulator_cm)-regulator_cm)
    return (math.hypot(masses[0],np.linalg.norm(state[6:9])) +
            math.hypot(masses[1],np.linalg.norm(state[9:12])) + string)


def propagate(parameters):
    """Run at most two pion decays; censor any leptonic decay (not yet implemented).

    No-recoil mode is deliberately an unphysical counterfactual.
    """
    p = parameters
    rng = np.random.default_rng(p.seed)
    b = benchmark(p.representation,p.mass_neutral_gev)
    m0 = p.mass_neutral_gev
    charged_mass = b["charged_mass_gev"]
    masses = [charged_mass,charged_mass]
    sigma = (p.lambda_ev*1e-9)**2
    p0 = charged_mass*p.beta/math.sqrt(1-p.beta**2)
    state = np.zeros(14)
    state[:3] = [p.impact_cm/2,0,5e-8]
    state[3:6] = [-p.impact_cm/2,0,-5e-8]
    state[6:9] = [0,0,p0]
    state[9:12] = [0,0,-p0]
    # Proper-time decay thresholds, independent for two charged endpoints.
    thresholds = (list(rng.exponential(b["lifetime_ps"],2))
                  if p.mode!="stable" else [math.inf,math.inf])
    charged = [True,True]
    events = []
    trace = []
    energies = []
    t = 0.0
    while t < p.duration_ps-1e-12:
        def rhs(_t, y):
            out = np.zeros(14)
            sep = y[:3]-y[3:6]
            f = -sigma/HBAR_GEV_PS * sep/math.hypot(np.linalg.norm(sep),p.regulator_cm)
            out[6:9] = f
            out[9:12] = -f
            for i in range(2):
                mom = y[6+3*i:9+3*i]
                en = math.hypot(masses[i],np.linalg.norm(mom))
                out[3*i:3*i+3] = C_CM_PER_PS*mom/en
                out[12+i] = masses[i]/en
            return out

        def make_event(i):
            def event(_t,y):
                return y[12+i]-thresholds[i]
            event.direction = 1
            event.terminal = True
            return event

        active = [i for i in range(2) if charged[i] and math.isfinite(thresholds[i])]
        energies.append(total_energy(state,masses,sigma,p.regulator_cm))
        sol = solve_ivp(rhs,(t,p.duration_ps),state,events=[make_event(i) for i in active],
                        method="DOP853",atol=1e-10,rtol=2e-9,
                        max_step=min(2.0,p.duration_ps/150))
        if not sol.success:
            raise RuntimeError(sol.message)
        trace.extend([float(ti),float(np.linalg.norm(yi[:3]-yi[3:6]))]
                     for ti,yi in zip(sol.t,sol.y.T))
        t = float(sol.t[-1])
        state = sol.y[:,-1].copy()
        energies.append(total_energy(state,masses,sigma,p.regulator_cm))
        choices = [j for j, ev in enumerate(sol.t_events) if len(ev)]
        if not choices:
            break
        i = active[choices[0]]
        if rng.random() >= b["pion_branching"]:
            events.append({"endpoint":i,"time_ps":t,"channel":"leptonic_unresolved"})
            break
        sl = slice(6+3*i,9+3*i)
        daughter, pion = pion_two_body(state[sl],masses[i],m0,rng)
        if p.mode=="recoil":
            state[sl] = daughter
        masses[i] = m0
        charged[i] = False
        thresholds[i] = math.inf
        events.append({"endpoint":i,"time_ps":t,"channel":"pion",
                       "pion_fourvector_gev":[float(pion[0]),*map(float,pion[1])],
                       "recoil_applied":p.mode=="recoil"})
    return {"model":p.representation,"mode":p.mode,"mass_neutral_gev":m0,
            "lambda_ev":p.lambda_ev,"lifetime_ps":b["lifetime_ps"],
            "delta_mev":b["delta_mev"],"pion_branching":b["pion_branching"],
            "events":events,"end_time_ps":t,"charged":charged,
            "segment_energies_gev":energies,
            "minimum_separation_including_production_cm":min(x[1] for x in trace),
            "maximum_separation_cm":max(x[1] for x in trace),
            "trace_time_ps_separation_cm":trace}
