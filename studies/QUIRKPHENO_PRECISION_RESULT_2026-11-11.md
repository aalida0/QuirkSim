# QuirkPheno precision confirmation — 2026-11-11

Verdict: PASS_SELECTED_STABLE_TRAJECTORY_RESOLUTION_COMPARISON.
The level-7/8 endpoint difference meets the selected 0.1-reference-length
comparison criterion over the sampled 0–130 ps interval. This supports starting
the weak-decay adapter; it is not a rigorous global error bound, validation of
full string invariants, or a recoil/annihilation result.

Input SHA-256: cabd22e0d82188a47931643405ed82844a2564c0e747eee054148e72573c40dc
Source: d91cb27ed3786bc63900faad892513841b5b6ed4.
Stable endpoint mass=225 GeV, Lambda=774 eV, point production, beta=.4 in
pair CM and beta_CM_z=.25. Not the split charged-parent benchmark or the
stretched-string conformal initial state.

## Completion and radial check

All eight executions and the build exit zero without timeout; macro hashes
match. One event, both endpoints, finite strictly timed steps, matching row
counts, kTime termination and at least 130 ps coverage in every case.
No G4Exception marker found in case logs.
Every endpoint has a px sign-change bracket including the radial turnaround
prediction near 111.4298 ps. At B=0 the level-8 endpoint brackets are
111.4224806–111.4324539 and 111.4258433–111.4329338 ps.
These bracket observations do not determine a precise magnetic root.

The B=0 analytic position residual is at most 2.795674e-9 mm, approximately
0.01097 reference lengths. Momentum residual norm at most 3.446e-5 MeV.
The plateau with refinement cautions against equating adjacent agreement
with absolute error. It is small for this selected criterion.
No total endpoint-plus-string energy-momentum history is available in ROOT.

## Magnetic precision

Worst endpoint difference at 2 T, cubic Hermite comparison on 2001 common times:

| Adjacent levels | Difference (nm) | Reference lengths |
|---|---|---|
| 5–6 | 0.0878787 | 0.344698 |
| 6–7 | 0.0439578 | 0.172421 |
| 7–8 | 0.0219692 | 0.0861724 |

Differences continue to halve. Both maxBoost and maxMerge change together,
so the scan does not identify a unique source of error.

Level 8 y displacement at 130 ps is -440.129033 and -440.129031 nm.
The level-7/8 endpoint-0 y change is 0.0214032 nm, about 0.00486% of
the transverse displacement. Macroscopic fractional accuracy and microscopic
encounter accuracy remain distinct questions.

Independent numerical sensitivity checks:
- 13001 common times reproduce max level-7/8 differences of
  2.1969187002e-8 and 2.1966337591e-8 mm.
- On that grid, Hermite vs position-only cubic spline differs by at most
  3.830e-12 mm at level 7 and 1.997e-12 mm at level 8.
These checks strongly reduce interpolation ambiguity for this comparison
but are not rigorous interpolation or solver error bounds.

Event CPU 0.024–0.186 s, excluding build/initialization/staging.
No extrapolation to a production ensemble or comparison with previous
long-time decaying trajectories is made.

## Frozen stable reference and next step

Use level 8 (maxBoost=0.0000390625, maxMerge=3.90625e-11 mm) as the selected
public-backend stable baseline, with level 7 as its nearby resolution control.
No additional stable refinement batch is required by these data for this
130 ps benchmark. Crossing beyond the approximately 222.86 ps radial return
has not been tested.

Proceed to the explicit parent/daughter endpoint handoff described in
QUIRKPHENO_PRECISION_AND_RECOIL_HANDOFF_2026-11-11.md.
That implementation is not present yet. It must preserve both histories and
partner ownership, change the endpoint species/mass/charge correctly, and
pass causal disturbance and local four-momentum checks. Adding a generic
decay table or resetting the string is insufficient.

No weak recoil result, physical survival probability or annihilation
suppression can be inferred. Doublet remains paused.
No further AF job was submitted during this audit.

Reproduce the main audit:
    python scripts/audit_quirkpheno_transport.py quirkpheno_precision.tar.gz --output audit.json

Derived data:
- studies/data/quirkpheno_precision_audit_2026-11-11.json
- studies/data/quirkpheno_precision_interpolation_2026-11-11.json
