# Variational string / weak recoil / radiation-envelope validation

**Project:** QuirkSim, independent open-source phenomenology implementation.
**Date:** 2026-11-10.
**Verdict:** `FINITE_ELEMENT_PROTOTYPE_PASS; CONTINUUM_WORLD_SHEET_NOT_VALIDATED; RADIATION_ANNIHILATION_PREDICTIONS_BLOCKED`.

This document records what was really implemented, numerically checked and
not accomplished. It must not be used to claim physical annihilation
suppression or a detector sensitivity projection.

## 1. Action-derived spatial string dynamics

The old two-particle central-force and Ng endpoint-force systems have no
independent string degrees of freedom. The new
[`quirksim/variational_string.py`](../quirksim/variational_string.py)
discretizes the classical Nambu-Goto action in laboratory-time gauge
into a piecewise-linear spatial embedding `x_j(t)`. Using midpoint
quadrature, its Lagrangian is

```text
L = -sum(endpoints) m_j sqrt(1 - |v_j|^2)
    -sum(interior) mu_g sqrt(1 - |v_j|^2)
    -T sum(segments i) ell_i sqrt(1 - |v_perp,i|^2)

T = Lambda_IC^2 / (hbar*c), ell_i = |x_{i+1}-x_i|,
v_mid,i = (v_i+v_{i+1})/2,
v_perp,i = v_mid,i - (v_mid,i . s_i) s_i.
```

The interior `mu_g` is a small **unphysical tangent-gauge regulator**,
and `ell_i` has a tiny explicit spatial regularization. Therefore this
is a spatial finite-element discretization, **not yet a certified continuum
Nambu-Goto worldsheet solution**. The published simple endpoint-force
equations of Li et al. 2020 are recovered only in appropriate local and
mesh-refined limits, not guaranteed at small vertex count.

For the semidiscrete action, the code builds the canonical mass matrix
`M=dP/dv` and convective term `dP/dx . xdot` and solves

```text
M dv/dt = c*dL/dx + q*v x B - (dP/dx)*xdot.
```

Unlike the previous approximation, string segments carry explicit
canonical **energy, linear momentum and angular momentum**. For a uniform
solenoid, mechanical angular momentum is not conserved because of the
external magnetic torque; without the external field, all three
Noether quantities must be conserved.

## 2. Verified short-time conservation and weak emission

Run: [finite-element mesh/weak-impulse gates](https://github.com/aalida0/QuirkSim/actions/runs/42838937110).

Six cases with **3, 5, 7 vertices** and interior gauge regulator fractions
`1e-2` or `1e-3`, initially length 0.1 cm and nonradial excitation.
Each case evolves a short stretch, then applies a single physical pion
decay through the finite-element **canonical** momentum, then continues.

- Before emission: absolute energy drift at most `1.71e-13 GeV`;
  canonical total momentum drift at most `2.85e-14 GeV`.
- At the pion impulse: absolute combined endpoint+string energy mismatch
  `4.43e-10 ... 3.98e-9 GeV`; momentum residual below
  `2.68e-11 GeV`; emitted pion energy near `0.172 GeV`.
- After emission: additional free-evolution energy drift
  below `1.14e-13 GeV`.
- The sum of string+endpoint angular momentum and emitted-pion
  orbital angular momentum is conserved numerically up to tiny absolute
  residuals relative to the huge macroscopic angular momentum.
- All six satisfy an explicit `1e-4 GeV` weak-energy gate; the
  production default is **fail closed** if the gate is violated.

These are **short-time numerical tests**, not a verified large-amplitude
multi-oscillation evolution. Fixed vertex-count and artificial gauge
mass can bias transverse wave propagation and self-intersection.

## 3. Transverse wave propagation and spatial refinement

Run: [linearized lowest transverse NG mode](https://github.com/aalida0/QuirkSim/actions/runs/42839273256).

Small-amplitude fundamental transverse excitation on a stretched string,
length 0.05 cm, with heavy but finite endpoints. The ideal fixed-end
linearized-NG half period is `L/c = 1.66782048 ps`; at this time the
string midpoint should have amplitude `-A`. The measured ratios:

| Nodes | midpoint transverse amplitude after L/c, divided by A | absolute deviation from -1 |
|---:|---:|---:|
| 5 | -0.979813 | 2.0187% |
| 7 | -0.992444 | 0.7556% |
| 11 | -0.995173 | 0.4827% |

The trend is consistent with a converging low-frequency transverse mode.
It does **not** demonstrate exact light-cone microcausality at finite
resolution, absence of numerical dispersion in higher modes, or global
continuum convergence in a folded moving string.

## 4. Actual event-driven proper-time weak decay

Run: [short triplet clocked pilot](https://github.com/aalida0/QuirkSim/actions/runs/42839828293).

Configuration: SU(2)L triplet, m_neutral=225 GeV, Lambda_IC=774 eV,
5 vertices, regulated pre-stretched length **0.1 cm**, beta_rel=0.4,
approximate pair x boost 0.25, uniform field 2T, observation interval
**5 ps**, chosen seed **25** to probe a rare early decay.

Sampled proper-time clocks are 12.07134 ps and 0.2021484 ps.
Second endpoint emits a pion at **lab t=0.25051456 ps**,
pion energy **0.15775513 GeV**. The finite-element total energy
mismatch across the canonical impulse is **8.46e-10 GeV**; the
momentum residual norm is ~4.29e-12 GeV. It continues dynamically
to 5 ps without violating the numerical action conservation gates.
The other endpoint remains charged during this observation window.

**This is a deliberately selected early-decay event, not a typical
triplet lifetime or a realistic LHC production sample.** It has
no post-decay close encounters before 5 ps, so none of the hazard
models is informed by that trajectory. The pre-stretched worldsheet
and approximate lab-time boost are not an exact pointlike hard
production initial-value prescription.

## 5. Radiation and annihilation: explicit assumption envelope ONLY

[`quirksim/encounter_hazards.py`](../quirksim/encounter_hazards.py)
is **not** self-consistent radiation dynamics. It provides
36 deliberately wide scenario combinations:

- Radius proxy `1/m_Q`, `1/sqrt(m_Q Lambda)`, or `1/Lambda`.
- Assumed independent per-geometric-contact annihilation probability
  `0, 0.01, 0.1, 1`.
- Assumed conditional per-contact glueball emission probability
  `0, 0.5, 1` with energy `k Lambda`.

It calculates conditional survival across contacts **before the first
weak decay**, to prevent automatically ignoring potentially annihilated
pairs. It also calculates a bookkeeping-only expected radiation energy;
**the string trajectory is never backreacted**. Independent contact
hazards and geometric radii are not calculated annihilation matrix
elements. The physically correct annihilation question involves
quantum wavefunction overlap, color/spin, and nonperturbative
deexcitation; a classical `r<1/Lambda` cut is not enough.

Li et al., *Detection prospects of long-lived quirk pairs at the
LHC far detectors*, Phys Rev D 109, 095005 (2024),
https://doi.org/10.1103/PhysRevD.109.095005,
use a phenomenological glueball emission probability epsilon per
oscillation with assumed energy scale Lambda. Our
per-**contact** diagnostic is **not identical** to that published
prescription: we require geometric proximity and never claim a
quantitatively calibrated emission rate. The nonperturbative
normalization is unknown.

## 6. Remaining blockers

**Cannot yet certify the publishable recoil–annihilation claim:**

1. Long-time evolved **continuum** worldsheet: mesh refinement and gauge
   regulator →0 under large stretch, folds, collisions, knots and boundary
   weak impulses; no physical remeshing algorithm currently.
2. Causal propagation and proper variational massive endpoint boundary
   behavior through repeated string crossings.
3. Self-consistent stochastic glueball and photon radiation, with
   conservation of emitted four-momentum and feedback on the worldsheet.
4. Quantum annihilation rates with realistic angular momentum, hidden
   group, spin, overlap and survival **prior** to the weak decay.
5. Correct pointlike pp-production initial conditions, realistic EW
   triplet channels and ISR kinematics from independent UFO/MG5.
6. Statistical uncertainties and convergence for ensembles.
7. Presently chosen early-decay seed and sub-first-return horizon
   do not show how recoil changes physical late-time observables.

**Verdict: numerical implementation substantial progress; paper
remains CONDITIONAL, NOT a validated GO.**

## Reproduce

```bash
git switch research/variational-worldsheet
python -m pip install -e '.[test]'
python -m pytest -q
python scripts/variational_gate_report.py
python scripts/validate_transverse_wave.py
python -m quirksim.worldsheet_experiment --seed 25 --duration-ps 5
```
