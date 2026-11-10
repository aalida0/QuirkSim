# Nambu–Goto constrained-gauge total-arclength refinement

**Date:** 2026-11-10  
**Model:** Dirac EW triplet benchmark, neutral mass 225 GeV, Λ_IC=774 eV, 2 T uniform field, beta_rel=0.4, approximate beta_CM,x=0.25, initial stretched-string length 0.1 cm. **Undecayed charged endpoints only.**

## Purpose and implementation

The previous long-time finite-element study stopped at an element-wise
fixed distance, causing the N=11 mesh to terminate around 227 ps while
N=5 and 7 reached ~459 ps. This was a **mesh-resolution-dependent stopping
criterion**, not evidence that the N=11 string physics differed.

The `research/variational-worldsheet` branch now computes the explicitly
mesh-independent **total spatial string arclength**

```text
L(t) = sum_{i=0}^{N-2} |x_{i+1}(t)-x_i(t)|
```

and stops when L(t)=L_cut, independently of vertex count or clustering.
This is a numerical safety boundary, not an annihilation radius or a
worldsheet reconnection/continuation prescription. Separate fail-closed
guards stop earlier if an individual element becomes underresolved or the
mesh length ratio becomes highly distorted.

Code:
- `quirksim/mesh_guard.py`: independent total-arclength event and
  separate mesh-quality events.
- `quirksim/equidistributed_ng.py`: strict fail-closed NG gauge solver.
- `scripts/equidistributed_gauge_stress.py`: explicit first-hitting-time
  diagnostic with physical-arclength residual.
- `scripts/check_total_arclength_convergence.py`: automated numerical
  convergence gates across 5/7/11 vertices at common cutoffs.
- `tests/test_total_arclength_gate.py`: tests including irregular
  interior gauge labels.

The solver assigns **zero artificial masses to internal string nodes**;
equal-arclength gauge constraints are solved simultaneously with the
semidiscrete Euler–Lagrange equations. This is still only a **midpoint
finite-element NG action approximation**, not a fully certified continuum
worldsheet.

## Results from 1,000 ps horizon requests

Every row below represents THREE independent resolutions: 5, 7 and 11
string vertices. They were evolved to the first crossing of their
common total-length threshold, **not for all 1,000 ps**.

| Common total L_cut (cm) | First threshold-crossing time (ps) | Max N=5,7,11 time spread (ps) | Max absolute energy drift (GeV) | Status |
|---:|---:|---:|---:|---|
| 1.0e-4 | ~227.208051 | 2.421e-6 | 3.41e-13 | Numerical gate PASS |
| 4.0e-5 | ~458.610359 | 2.261e-7 | 3.98e-13 | Numerical gate PASS |
| 1.0e-5 | ~458.611638 | 3.016e-7 | 3.98e-13 | Numerical gate PASS |
| 2.0e-6 | ~458.611979 | 1.199e-6 | 3.98e-13 | Numerical gate PASS |

At each cutoff the total-arclength event residual is below
`3e-15 cm` and the equal-arclength gauge remains exceptionally
well balanced up to the stop. A different threshold can select a
**different close encounter**: the larger 1e-4 cm cutoff stops at
~227 ps; the 4e-5 cm and smaller cutoffs reach a later encounter at
~459 ps. These are threshold crossings, not proof of a mathematically
zero-length singularity.

Public GitHub Actions:
- [2 common cutoff + unit-test physics gate](https://github.com/aalida0/QuirkSim/actions/runs/42850022362)
- [2 tighter total-cutoff convergence checks](https://github.com/aalida0/QuirkSim/actions/runs/42850649116)
- [independent radial massive-endpoint exact 1D benchmark](https://github.com/aalida0/QuirkSim/actions/runs/42844612951)

## Crucial unsolved physical-resolution question

At the tightest tested **total string arclength L_cut=2e-6 cm**, the
endpoint chord lengths are

| Vertices | Endpoint separation (cm) |
|---:|---:|
| 5 | 1.998061959e-6 |
| 7 | 1.993005809e-6 |
| 11 | 1.970712308e-6 |

The 5→11 endpoint separation difference is about
**2.73e-8 cm**. The reference confinement length
`hbar*c/Lambda_IC ≈ 2.55e-8 cm` is of the SAME ORDER. Thus the
**stopping time and total arclength converge**, but the transverse
geometry at an extremely close encounter is **NOT demonstrated to
converge below the confinement scale**. The code now reports a
separate core-reference geometry gate, marked **BLOCKED** if the
inter-resolution chord spread exceeds 10% of this reference length.

Neither `L_cut` nor `1/Lambda` is itself a calculated probability
for quirk annihilation. Moreover, this study uses a pre-stretched
string slice and a fixed boost, not realistic pointlike pp+jet
production, and does not apply weak recoils to the NEW gauge solver.

## Current physics verdict

- **PASS:** Explicit total-arclength cutoff independent of mesh
  discretization; near-identical pre-stop times and global energy
  conservation at 5, 7, 11 vertices; no artificial internal node masses.
- **BLOCKED:** Core-scale geometrical convergence of the endpoint
  separation near very close returns; existence/uniqueness of
  physically correct worldsheet continuation or reconnection through
  zero length; finite causal weak recoil on the NEW gauge; radiation,
  quantum annihilation, and survival before weak decay.
- **DO NOT CLAIM:** Completed 1,000 ps physical evolution, annihilation
  suppression, collider reach, or a validated continuum 3D
  Nambu–Goto result.

Next numerical target: mesh and timestep refinement of the *transverse
chord* and field/string momentum exchange at the closest approach,
then a separately motivated continuation prescription. Increasing
the number of Monte Carlo decay samples before these gates pass is
not useful.

### Reproduction

```bash
git switch research/variational-worldsheet
python -m pip install -e '.[test]'
python -m pytest -q tests/test_total_arclength_gate.py
python scripts/check_total_arclength_convergence.py
python scripts/check_total_arclength_convergence.py --fractions 0.0001,0.00002
```
