# Specified vector-current production and local weak-decay spin maps

## New result

For a Dirac charged pair produced by a pure heavy vector current in
q qbar -> gamma*/Z* -> Q+ antiQ+, the canonical spin density in the
production pair centre-of-mass frame has **zero spin-singlet occupation**.
This holds for either transverse source helicity separately and for any
nonnegative helicity mixture. The final heavy Z current is vectorlike for
the Y=0 vectorlike triplet; an initial-quark axial coupling changes source
helicity weights, not this conclusion. No massive-initial-quark, extra
hard-jet, loop, or axial heavy-current amplitudes are included here.

This is a specified production-model calculation, not information recovered
from the AF gun: the gun generator has not implemented a production spin
density. The existing compact spin-zero widths must not be weighted by
an assumed 1/4 singlet population for this vector-current source.

## Derivation and independent checks

Use canonical physical spins for the outgoing particle and antiparticle,
with spinor normalizations ubar u=2m, vbar v=-2m. In the antiparticle
spinors eta_s=(-i sigma_2) chi_s*, so both physical spin indices transform
in the same SU(2) convention. The four-spin amplitudes are evaluated as
J_i(s,t)=ubar(p,s) gamma_i v(-p,t). Contracting J with either circular
transverse source gives a vector orthogonal to
(|up down>-|down up>)/sqrt(2). Equivalently, the vector source produces
the spin-triplet partial waves (3S1 and 3D1); it does not produce 1S0.

Tests cover beta=0,0.4,0.9, three outgoing directions and each source
helicity. The spin-summed amplitudes reproduce the independent Born angular
dependence 1+cos^2(theta)+(1-beta^2)sin^2(theta). This is an algebraic and
numerical derivation, not an externally reproduced event-generator result.
Attempts to obtain additional primary literature during this calculation
failed with search-service 503 errors; no new literature confirmation is
claimed for this step.

## Pion decay does not randomly reset the spin

The leading vector-current pion amplitude has
q_mu ubar_d gamma^mu u_parent = (M-m) ubar_d u_parent by the on-shell Dirac
equations. In the parent rest frame, this scalar bilinear is proportional
to the identity in canonical physical spin labels, independently of the
pion direction. In a general frame its normalized 2x2 matrix is unitary:
it gives a canonical spin rotation, not depolarization. The physical
antiparticle map agrees up to an irrelevant overall phase.

The implementation computes both the contracted vector current and scalar
reduction independently on sampled decay kinematics. It also checks the
antiparticle map and inverse Lorentz-boost rotations. All four focused
tests pass. The eight logged pion decays in the six-case 260 ps archive
give local canonical rotation angles at most 5.0641e-5 radians. These
local maps cannot simply be applied to the production density across a
missing propagation history.

At the fixed beta_pair=0.25 boost, the production density's **lab canonical
spin-only** singlet fraction is 0.0015256284. Transforming back to the
production pair CM restores zero. This frame-dependent spin-only number
is not an annihilation occupation: spin and momentum must be transformed
together into a consistently defined matching frame. This check prevents
an artificial spin-singlet population from a basis change alone.

## Consequence for the annihilation task

The previous spin-zero compact calculation remains a valid specified-sector
subcalculation. Its occupation is not established for the actual source.
For a strictly spin-conserving nonrelativistic propagation model with this
DY source, that spin-zero contribution would remain absent. This does not
imply that the pair is stable or that the full annihilation rate vanishes:
spin-one, higher partial waves, spin-changing interactions and radiative
transitions must be considered.

Continuous spin transport under string and electromagnetic forces has not
been supplied by the spinless endpoint histories. Local decay maps alone
do not close it. A spatial production density and radiation/capture model
are also still required to predict with/without-recoil probabilities.
No physical total annihilation probability is assigned by this result.

## Finite-velocity neutral hard scattering

The new `neutral_hard_scattering.py` also evaluates the full vectorlike
t/u amplitudes for a bare neutral Dirac particle/antiparticle pair at finite
CM velocity, retaining all four initial physical spin states and final vector
polarizations. The WW neutral s-channel vanishes because T3=Y=0. Fermion
masses are degenerate at this tree order; the IC-singlet color projections
are N_IC for WW and (N_IC^2-1)/(4 N_IC) relative to the abelian two-photon
kinematic kernel at the same gauge coupling. The latter includes the
identical-final-boson factor.

For m0=225 GeV, N_IC=2, alpha_IC,hard=0.04 and beta_each=0.4:

| Channel | Pure spin-singlet sigma*v (GeV^-2) | Isotropic spin-triplet average sigma*v (GeV^-2) |
| --- | ---: | ---: |
| WW | 8.53644e-7 | 1.25351e-7 |
| g' g' | 1.18203e-7 | 1.52387e-8 |

Here v is the CM Moller flux velocity 2 beta, not the invariant relative
speed 2 beta/(1+beta^2). At beta=0 the singlet limits reproduce the
previous analytic threshold kernels and the triplet terms vanish. The
finite-velocity spin-triplet terms include higher orbital partial waves.
They must **not** be used as l=0 bound-state widths by multiplying them
by an S-wave contact density.

Angular integration changes by less than 1e-8 relative between the two
checked quadratures at beta=0.4. The massless-vector spin average also
reproduces the analytic inverse Breit-Wheeler QED result at beta=0.05,
0.4 and 0.8. Three additional tests check this reference, threshold
projection, positive spin kernels and angular convergence. The code
outputs the full 4x4 hard spin matrices so a supplied incoming density
can be used in a consistently defined spin and momentum frame.

These hard coefficients resolve finite-velocity bare two-body spin coverage
for these two channels. They do not supply the incoming spin/orbital
wavefunction, rescattering, higher-multiplicity channels or the quantum
population/transition data needed for an annihilation probability.

Reproduce:

```bash
python scripts/calculate_neutral_finite_velocity.py --output neutral_finite_velocity.json
PYTHONPATH=. python -m unittest discover -s research/quirkpheno_recoil -p test_neutral_hard_scattering.py
```

Reproduce the local-log audit and source-density calculation:

```bash
python scripts/audit_triplet_production_spin.py triplet_first_return_audit.tar.gz \
  --output triplet_production_spin.json
PYTHONPATH=. python -m unittest discover -s research/quirkpheno_recoil -p test_production_spin.py
```

The output includes the original archive SHA-256, production density and
all eight local rotation matrices. No new AF run was used.
