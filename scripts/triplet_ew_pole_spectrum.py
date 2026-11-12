#!/usr/bin/env python3
"""One-loop electroweak pole-mass *difference* for an SU(2)L Y=0 fermion triplet.

Minimal triplet with one common gauge-invariant tree mass M, no Yukawa
mixing or additional non-SM mass corrections.

For charges +/-1 relative to neutral 0, the finite, gauge-invariant
one-loop difference is [Cirelli, Fornengo, Strumia, hep-ph/0512090]
    DeltaM = alpha2 * M/(4*pi) [f(mW/M)-cW^2*f(mZ/M)-sW^2*f(0)]
    f(r) = 2 * int_0^1 (1+x) log[x^2+(1-x)r^2] dx, f(0)=-5.
This holds for both vectorlike Dirac and minimal Majorana triplet
gauge interactions; it does NOT calculate a complete EW-renormalized
collider amplitude, and its inputs are renormalization-scheme dependent
at higher order.

Default numerical inputs are the SM.fr / MG5-style G_F,M_Z,alpha(M_Z)
tree relations: MW is *derived* from the other three. This scheme gives
a mass-dependent result and the ~166 MeV asymptotic limit. Using
physical MW and an independently specified G_F corresponds to a
different finite-order scheme. Record the inputs when quoting results.

Pole-mass convention: define the neutral *pole* mass as the reference M,
then mcharged(pole) = M + DeltaM. This definition is NOT meant to
insert a gauge-breaking tree mass difference into the NLO-QCD UFO.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

PI = math.pi


def _simpson(func, a: float, b: float) -> float:
    c = (a+b)/2
    return (b-a)/6*(func(a)+4*func(c)+func(b))


def _adaptive_simpson(func, a: float, b: float, whole: float,
                      eps: float, depth: int) -> float:
    mid=(a+b)/2
    left=_simpson(func,a,mid)
    right=_simpson(func,mid,b)
    delta=left+right-whole
    if depth <= 0:
        if abs(delta)> 16 * eps:
            raise ArithmeticError("EW loop integral failed precision convergence")
        return left+right+delta/15
    if abs(delta) <= 15*eps:
        return left+right+delta/15
    return (_adaptive_simpson(func,a,mid,left,eps/2,depth-1)
            +_adaptive_simpson(func,mid,b,right,eps/2,depth-1))


def loop_f(r: float) -> float:
    """Finite one-loop function, integral evaluated without scipy/numpy."""
    if not math.isfinite(r) or r < 0:
        raise ValueError("Loop argument r must be finite and nonnegative")
    if r == 0:
        return -5.0  # exact integral
    # x=t^2 removes the endpoint structure from log(x^2) as r->0.
    def func(t):
        x=t*t
        return 4*t*(1+x)*math.log(x*x+(1-x)*r*r)
    return _adaptive_simpson(func,0.0,1.0,_simpson(func,0.0,1.0),3e-11,24)


@dataclass(frozen=True)
class Spectrum:
    input_mass_gev: float
    neutral_pole_mass_gev: float
    charged_pole_mass_gev: float
    delta_m_gev: float
    delta_m_mev: float
    mw_gev: float
    mz_gev: float
    alpha_em: float
    alpha2: float
    sin2theta_w: float
    scheme: str


def triplet_spectrum(
    m_gev: float,
    *,
    aewm1: float=127.9,
    gf: float=1.1663787e-5,
    mz_gev: float=91.1876,
    mw_gev: float | None = None,
) -> Spectrum:
    vals=(m_gev,aewm1,gf,mz_gev)
    if not all(math.isfinite(x) and x>0 for x in vals):
        raise ValueError("All input mass, coupling and G_F values must be positive finite")
    alpha=1./aewm1
    if mw_gev is None:
        # Same tree-level G_F-MZ-alpha input scheme as the standard SM.fr
        discriminant=mz_gev**4/4 - PI/math.sqrt(2)*alpha/gf*mz_gev*mz_gev
        if discriminant<0:
            raise ValueError("Inconsistent G_F,MZ,alpha: complex tree-level MW")
        mw=math.sqrt(mz_gev*mz_gev/2+math.sqrt(discriminant))
        scheme="GF_MZ_ALPHA_MZ_TREE"
    else:
        if not math.isfinite(mw_gev) or mw_gev<=0:
            raise ValueError("Input MW must be positive finite")
        mw=mw_gev
        scheme="MW_MZ_ALPHA_INPUT"
    sw2=1-(mw/mz_gev)**2
    if not 0<sw2<1:
        raise ValueError("Physical electroweak weak angle requires 0 < MW < MZ")
    cw2=1-sw2
    alpha2=alpha/sw2
    combo=loop_f(mw/m_gev)-cw2*loop_f(mz_gev/m_gev)-sw2*loop_f(0)
    dm=alpha2*m_gev/(4*PI)*combo
    if dm<=0 or not math.isfinite(dm):
        raise ValueError("Charged-neutral splitting is nonpositive or nonfinite")
    return Spectrum(
        input_mass_gev=m_gev,
        neutral_pole_mass_gev=m_gev,
        charged_pole_mass_gev=m_gev+dm,
        delta_m_gev=dm,
        delta_m_mev=1000*dm,
        mw_gev=mw,
        mz_gev=mz_gev,
        alpha_em=alpha,
        alpha2=alpha2,
        sin2theta_w=sw2,
        scheme=scheme,
    )


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mass",type=float,default=225.,help="Neutral pole mass, GeV")
    p.add_argument("--aewm1",type=float,default=127.9,help="Inverse electromagnetic coupling")
    p.add_argument("--gf",type=float,default=1.1663787e-5,help="Fermi constant GeV^-2")
    p.add_argument("--mz",type=float,default=91.1876,help="Z pole mass in GeV")
    p.add_argument("--mw",type=float,help="Optional W mass; changes scheme if supplied")
    p.add_argument("--output",type=Path,help="Write JSON spectrum to this path")
    a=p.parse_args()
    result=triplet_spectrum(a.mass,aewm1=a.aewm1,gf=a.gf,
                           mz_gev=a.mz,mw_gev=a.mw)
    print(f"Neutral pole mass: {result.neutral_pole_mass_gev:.9f} GeV")
    print(f"Charged pole mass: {result.charged_pole_mass_gev:.9f} GeV")
    print(f"DeltaM_EW one-loop: {result.delta_m_mev:.6f} MeV")
    print(f"Input scheme: {result.scheme}, MW={result.mw_gev:.6f} GeV, "
          f"sin2thetaW={result.sin2theta_w:.6f}")
    print("CAUTION: one-loop spectrum only; NOT a UV/R2-complete NLO-EW UFO.")
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(asdict(result),indent=2)+"\n")
        print("Saved:",a.output)


if __name__=="__main__":
    main()
