# Nambu–Goto near-contact spatial quadrature audit

**Date:** 2026-11-10.  
**Status:** MIDPOINT_APPROXIMATION_INADEQUATE; GAUSS2/GAUSS3_AGREEMENT_PASS;
**SELECTED_PRE_CONTACT_CORE_REFERENCE_GEOMETRY_PASS** (13/17/21 vertices);
FULL_CONTACT_CONTINUATION_AND_ANNIHILATION_BLOCKED.

## Motivation

In the research constrained-gauge solver, strings with 5,7,9,11,13,15,17
vertices have nearly identical first-hitting times of a fixed **total
arclength** L_cut=2.0e-6 cm, but their endpoint **chord length** is not
monotonic with vertex count. At n=9 vs n=5 the discrepancy is order
hbar*c/Lambda_IC. This cannot be fixed by smaller ODE steps or a smaller
short-distance regulator.

The previous action used one midpoint evaluation of the nonlinear
Nambu–Goto segment integrand: T * length *
sqrt(1 - v_perp(midpoint)^2). Although x(s) and v(s) are linear within
a segment, sqrt(1-v_perp(s)^2) is **not linear**. One midpoint evaluation
is thus an uncontrolled spatial quadrature for large transverse
velocity gradients.

## Independent action-derived discretizations

`quirksim/ng_quadrature.py` now implements Gauss–Legendre integration
of the segment Nambu–Goto action with **1, 2 or 3 integration points
per element**. Canonical momenta P=dL/dv, mass matrix M=dP/dv,
position force dL/dx, convective contribution dP/dx * xdot, and
energy all use **the same quadrature**. This is not an ad-hoc force
correction: the semidiscrete equations are Euler–Lagrange equations
of an explicitly defined discrete action.

`tests/test_ng_quadrature.py` checks midpoint equivalence, independent
finite-difference gradients of the action, and short-time magnetic
energy conservation. The GitHub CI tests passed.

## Same near-contact slice

Same point in parameter space as earlier: triplet m_neutral=225 GeV,
Lambda_IC=774 eV, pre-stretched length 0.1 cm, Bz=2 T,
beta_rel=0.4, beta_CM,x=0.25; no weak decays or radiation.

All numbers are endpoint chord length in **cm** when the total
string arclength first reaches 2.0e-6 cm.

| Vertices | 1-point midpoint | 2-point Gauss | 3-point Gauss |
|---:|---:|---:|---:|
| 5 | 1.9980617e-6 | **1.9412613e-6** | 1.9414174e-6 |
| 9 | 1.9401365e-6 | **1.9758246e-6** | 1.9758368e-6 |
| 13 | 1.9893920e-6 | **1.9979658e-6** | 1.9979662e-6 |

**Two- and three-point Gauss agree far better than either agrees with
the midpoint approximation**. For n=5, the midpoint→Gauss2 chord shift
is ~5.68e-8 cm (~2.23 times hbar*c/Lambda_IC), whereas the
Gauss2→Gauss3 difference is ~1.56e-10 cm. For n=9, midpoint→Gauss2
is ~3.57e-8 cm and Gauss2→Gauss3 ~1.22e-11 cm. For n=13,
Gauss2→Gauss3 is ~4.42e-13 cm.

Consequently, **midpoint quadrature was a major cause of the
previous numerical uncertainty**. However, the improved Gauss2
result still changes by ~5.67e-8 cm between n=5 and n=13;
higher-order quadrature does not by itself establish spatial
continuum convergence.

Conserved energy differs by at most ~2e-12 GeV over these test
histories; energy conservation cannot replace an explicit
endpoint-geometry convergence test.

GitHub CI:
- [Gauss 2/3 benchmark, 5/9/13 vertices](https://github.com/aalida0/QuirkSim/actions/runs/42893666233)
- [Independent midpoint spatial-refinement anomaly](https://github.com/aalida0/QuirkSim/actions/runs/42892811180)

## Final independent 13/17/21-vertex validation: PASS for this selected slice

A **new independent CI workflow** runs all four physics configurations
in separate jobs, downloads the actual per-run JSON artifacts, and then
runs `scripts/validate_ng_gauss_refinement.py` as a fail-closed
aggregate geometry test.

Common inputs: Bz=2 T, pair boost=0.25, Lambda_IC=774 eV,
neutral mass=225 GeV, initial stretched length=0.1 cm,
Gauss–Legendre order=2, short-distance regulator
**2.0e-10 cm**, total arclength cutoff **2.0e-6 cm**.

| Mesh vertices | Gauss-2 endpoint chord (cm) |
|---:|---:|
| 13 | 1.997965607203e-6 |
| 17 | 1.999469084312e-6 |
| 21 | 1.999712177476e-6 |

- Largest mesh chord spread: **1.74657027e-9 cm**,
  **0.0685079 × (hbar*c/Lambda_IC)**.
- Required convergence diagnostic: **< 0.1 × reference length — PASS**.
- Gauss-2 vs Gauss-3 at 21 vertices: **1.77547e-14 cm**, PASS.
- 13/17/21 total-length hitting-time spread:
  **7.34e-8 ps**, PASS.
- All four configurations reached the shared **total-arclength cutoff**,
  with no per-element mesh stop.
- These checks deliberately **stop before the unresolved gauge pinch**.
  The requested 1000 ps was not physically completed.
- [Passing CI and four raw JSON trajectories' endpoint summaries](https://github.com/aalida0/QuirkSim/actions/runs/42895003435).

**Interpretation:** The higher-order action and sufficient spatial
resolution substantially mitigate the previous midpoint aliasing.
This establishes pre-contact geometric stability in one controlled
benchmark, *not* full continuum Nambu–Goto validation. A 10% tolerance
relative to a confinement reference length is a chosen **numerical
diagnostic**, not an annihilation error budget or a quantum calculation.

## Remaining physics gates

1. Extrapolate the Gauss2 (and a cross-check Gauss3) endpoint geometry
   toward finer spatial meshes. The core-reference criterion is
   |r_min(N)-r_min(N_ref)| < 0.1 hbar*c/Lambda_IC and should include
   regulator/timestep sensitivity; this is a diagnostic only.
2. Develop well-posed worldsheet continuation through near-contact
   gauge degeneracy and test discrete conservation and causal
   propagation, rather than manufacturing a crossing event.
3. Revalidate weak-decay endpoint impulse in the new zero-interior-mass
   Gauss-constrained action; old finite-element results do not transfer
   automatically.
4. Only after these: radiation and quantum annihilation envelopes with
   realistic pre-decay survival.

**No physical annihilation rate, detector simulation or LHC
sensitivity is implied by these results.**
