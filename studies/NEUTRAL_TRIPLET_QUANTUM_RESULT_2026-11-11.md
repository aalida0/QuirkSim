# Neutral triplet: coupled spin-zero annihilation and unresolved matching

## What is now calculated

A specified compact, neutral-connected 1S spin-zero eigenstate has been solved
with charged-neutral flavor mixing. This is a new model calculation, not
application of the charged benchmark's widths to neutral endpoints.

The field representations match the Dirac (1,3)_0, dark fundamental example
in Buen-Abad, Marques-Tavares and Schmaltz,
[PRD 92 (2015) 023531, Sec. II](https://doi.org/10.1103/PhysRevD.92.023531).
Its cosmological dark-radiation regime and thermal color-averaged rate are
not transferred to a confined IC-singlet quirkonium. The Majorana wino
calculation in [Asadi et al., JCAP 02 (2017) 005](https://doi.org/10.1088/1475-7516/2017/02/005)
is useful background on coupled-channel capture; its identical-particle
spin restrictions and normalization are not assumed for this Dirac pair.

### Specified model and potential

Use m0=225 GeV, delta=m_charged-m0=0.157170276 GeV, N_IC=2,
alpha_IC,potential=alpha_IC,hard=0.04 and sqrt(tension)=774 eV.
Couplings are explicit fixed LO inputs; no Lambda_MS/string-tension match is
inferred. Electroweak inputs use GF, MW, MZ consistently, with
alpha2=sqrt(2) GF MW^2/pi and alpha_EM=alpha2 sin^2(thetaW).

The Dirac Q=0 flavor basis is (+ anti+, 0 anti0, - anti-), without
identical-neutral normalization factors. The central potential is

V = (-C_F alpha_IC/r + tension*r) I +
[[2 delta-A, -B, 0], [-B, 0, -B], [0, -B, 2 delta-A]],

A=(alpha_EM+alpha2 cos^2(thetaW) exp(-MZ r))/r,
B=alpha2 exp(-MW r)/r.

This is derived from the triplet generators with a phase convention making
the off-diagonals negative. At short distance its EW eigenvalues are
alpha2*(-2,-1,+1)/r, as required for 3 x anti3 = 1+3+5.
Spin-dependent potentials, radiative capture, running and higher-order
potential corrections are omitted from this LO static model.

### Hard amplitude matrices

Define K_X by Gamma_X=psi(0)^dagger K_X psi(0). All five matrices are
positive rank-one matrices, retaining coherent flavor interference.
Set b=pi/m0^2, a=alpha2, e=alpha_EM, c=cos^2(thetaW),
F_V=(1-MV^2/m0^2)^(3/2)/(1-MV^2/(2m0^2))^2 above threshold, zero below.

| Final state | Coefficient multiplying v v^dagger | Flavor amplitude v |
| --- | --- | --- |
| g' g' | b alpha_IC,hard^2 (N_IC^2-1)/N_IC | (1,1,1) |
| W+ W- | 8 b N_IC a^2 F_W | (1/2,1,1/2) |
| gamma gamma | 4 b N_IC e^2 | (1,0,1) |
| gamma Z | 8 b N_IC e a c [1-MZ^2/(4m0^2)] | (1,0,1) |
| Z Z | 4 b N_IC (a c)^2 F_Z | (1,0,1) |

The gamma-Z expression is set to zero below its threshold. Tree amplitudes
use degenerate fermion masses; the loop splitting is retained in the
nonrelativistic channel thresholds. This avoids inserting an unmatched
symmetry-breaking mass into selected tree propagators. The spin-zero
vector current does not contribute to s-channel fermion production; no
Higgs Yukawa interaction is included in this vectorlike model.

Independent explicit Dirac spinors, t/u propagators and physical vector
polarizations reproduce the neutral WW and charged neutral-boson diagonal
coefficients at m=100,225,1000 GeV. Those threshold spin matrices have rank
one; the other three spin eigenvalues vanish to numerical roundoff. Flavor
amplitudes follow the diagonal triplet generator products. The IC factor
uses the previously checked singlet color trace and spin projection.
This validates the leading S-wave spin-zero two-body sector. It is not
a full finite-velocity or spin-one annihilation calculation.

### Compact eigenstate result

The coupled radial Schrodinger equation has kinetic term -u''/m0,
u(0)=u(rmax)=0 and integral sum_i |u_i|^2 dr=1. Meshes 4000, 8000,
16000, 32000 were checked. The contact amplitudes are extrapolated from
u_i/r at the first five radii. Final partial-width changes are below 0.05%.
The solver also reproduces analytic Coulomb energy and contact density
in an EW-off control to 1e-4 relative precision. These are numerical
checks, not estimates of omitted perturbative uncertainty.

| Channel | Width of this specified LO compact eigenstate (keV) |
| --- | ---: |
| W+ W- | 16.6050 |
| g' g' | 2.81546 |
| gamma gamma | 0.0176381 |
| gamma Z | 0.117831 |
| Z Z | 0.194076 |

hbar divided by this modeled two-boson sum is 3.33e-20 s. The
eigenstate occupation in the transported pair remains unspecified. These
numbers are not assigned to the first-return trajectories.

## Why matching and a recoil-dependent probability are not identified

The saved gun/weak-decay files provide classical sample means. They do not
specify a spin density matrix, a quantum wavepacket or a radiation history.
Even two individually unpolarized spins can be in a singlet, an isotropic
triplet mixture or a product mixture: the one-particle density matrices
are identical but their singlet occupations are 1, 0 and 1/4 respectively.
The unpolarized weak-decay sampler cannot distinguish these states.
For a specified vector-current Drell-Yan source the spin state must instead
be derived from the production amplitude and transported through decay;
it cannot be replaced by an arbitrary 1/4 singlet occupation. The gun has
not implemented that production amplitude.

There is an independent spatial ambiguity. In a *ballistic Gaussian contact
ansatz*, keeping identical classical position/momentum means but varying
the transverse width w gives the contact ratio
R=exp[-b^2/(2w^2)] for impact b relative to zero impact. At the observed
b=4.154395 micrometres, log(R) is about -8.63e6 for w=1 nm, -8.63 for
w=1 micrometre and -0.0863 for w=10 micrometres. This diagnostic is not a
physical uncertainty band: it deliberately changes an unspecified input,
and excludes focusing, capture and radiation. It demonstrates why the
classical miss alone does not determine a quantum annihilation probability.

Finally, a local radiation/capture kernel must determine energy, momentum,
spin and angular-momentum transfer. An average glueball efficiency alone
does not supply this quantum transition law. The effects cannot be covered
by a unique numerically justified radiation uncertainty band from the
available histories.

## Status of the three requested items

1. Leading spin-zero S-wave two-body annihilation is implemented and checked
   for a specified coupled compact eigenstate. Spin-one, finite-energy and
   higher-order coverage remains incomplete.
2. Unique quantum matching from the existing classical histories is not
   determined. A specified production density operator and an evolution law
   for weak decay, spin, radiation and capture are required.
3. A physical with/without-recoil annihilation probability is consequently
   still UNKNOWN. Neither the classical radius nor the compact-state width
   supplies the missing occupation and transition probabilities.

Release: the declared compact spin-zero model calculation can be used as a
subcalculation. The full physical recoil suppression claim is not released.
No extra AF transport run is needed to resolve the missing quantum inputs.

Reproduce:

```bash
python scripts/calculate_neutral_triplet.py --output neutral_triplet_quantum.json
PYTHONPATH=. python -m unittest discover -s research/quirkpheno_recoil -p test_neutral_triplet.py
```
