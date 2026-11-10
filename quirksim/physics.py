"""Leading-order electroweak doublet/triplet mass splittings and weak widths.

Masses/widths are in GeV; lifetimes in ps. No precision electroweak claims.
Thomas & Wells, arXiv:hep-ph/9804359; standard Y=0 triplet loop kernel.
"""
import math
from scipy.integrate import quad

MZ = 91.1876
MW = 80.377
GF = 1.1663787e-5
ALPHA_INV_MZ = 127.95
VUD = 0.97373
FPI = 0.1302
MPI = 0.13957039
ME = 0.00051099895
MMU = 0.1056583755
HBAR_GEV_PS = 6.582119569e-13


def doublet_splitting(m):
    """Thomas-Wells Eqs.(2-3), one-loop splitting in GeV."""
    if m <= 0:
        raise ValueError("mass must be positive")
    q = m / MZ
    upper = q / (1 + q)
    def kernel(y):
        if y == 0:
            return 0.0
        t = y / (1 - y)
        return ((1 + t / q) * math.log1p((1 - t / q) / (t*t))
                / ((1 - y)**2 * math.pi))
    integral = quad(kernel, 0, upper, epsabs=1e-10, epsrel=1e-10, limit=300)[0]
    return MZ * integral / (2 * ALPHA_INV_MZ)


def triplet_splitting(m):
    """Y=0 triplet electroweak one-loop splitting, on-shell GF scheme."""
    if m <= 0:
        raise ValueError("mass must be positive")
    s2 = 1 - (MW/MZ)**2
    alpha2 = math.sqrt(2) * GF * MW**2 / math.pi
    def f(z):
        return quad(lambda x: 2*(1+x)*math.log(x*x + (1-x)*z*z),
                    0, 1, epsabs=1e-10, epsrel=1e-10, limit=300)[0]
    return alpha2*m/(4*math.pi) * (f(MW/m) - (1-s2)*f(MZ/m) - s2*f(0))


def phase(b):
    if not 0 <= b < 1:
        raise ValueError("phase-space ratio outside [0,1)")
    if b < 1e-8:
        return 1.0
    a = math.sqrt(1-b*b)
    return 1-4.5*b*b-4*b**4 + 7.5*b**4/a*math.atanh(a)


def widths(delta, representation):
    """Leading pion and lepton partial widths; CG^2 triplet/doublet = 2."""
    if representation not in ("doublet", "triplet") or delta <= 0:
        raise ValueError("bad representation or delta")
    mult = 2 if representation == "triplet" else 1
    pion = 0.0
    if delta > MPI:
        pion = GF**2 * VUD**2 * FPI**2 * delta**3 / math.pi * math.sqrt(1-(MPI/delta)**2)
    def leptonic(ml):
        if delta <= ml:
            return 0.0
        b = ml/delta
        return GF**2 * delta**5 / (15*math.pi**3) * math.sqrt(1-b*b) * phase(b)
    return {"pion": mult*pion, "electron_neutrino": mult*leptonic(ME),
            "muon_neutrino": mult*leptonic(MMU)}


def benchmark(representation="triplet", mass_gev=225.0):
    delta = (doublet_splitting(mass_gev) if representation == "doublet"
             else triplet_splitting(mass_gev))
    w = widths(delta, representation)
    total = sum(w.values())
    return {"representation": representation, "neutral_mass_gev": mass_gev,
            "charged_mass_gev": mass_gev+delta, "delta_mev": 1000*delta,
            "widths_gev": w, "lifetime_ps": HBAR_GEV_PS/total,
            "pion_branching": w["pion"]/total}
