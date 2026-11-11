# Event-side consistency fix and rejected boundary-method pilots

## Narrow fix accepted

The characteristic event handling and `next_jump` traversal use a 32-ULP
roundoff window, whereas `UnitRayHistory.value` used 16 ULPs. At a coordinate
17–31 ULPs before a front, traversal therefore regards the front as processed
but lookup can return its pre-front direction. Reconstructing the outgoing
massive-boundary ray at that same instant disagrees with an already stored
post-event reflected ray, even for a zero-length trial.

The main research branch changes the lookup window to the same existing
32-ULP policy. It does not change force laws, acceptance tolerances, sampling
controls, decay probabilities or quantum matching. A regression reproduces
the old inconsistency at 24 ULPs, and checks 17/24/31 ULPs after the fix.
Outside the window, lookup remains pre-front; the exact left limit remains
available; a jump still has zero integrated impulse. The selected legacy
characteristics, impulse, matching and diagnostic tests pass: 26 tests.

This demonstrates and fixes a real side-policy inconsistency. It does **not**
prove that the reported 315 ps failure was caused by it. The 15 reflected
fronts inspected in the available old 39 ps checkpoint were within -2 to +4
ULPs of their event coordinates, not in the inconsistent window. No later
full-history checkpoint from the three-hour run is available. A fresh bounded
prefix with the fix reached 27.1076 ps and stopped at the 120-second diagnostic
budget; it does not validate the 315 ps case or establish a speed improvement.

## Alternatives tested and not promoted

Separate development pilots tested two larger changes:

1. An implicit conservative boundary update using common incoming/outgoing
   characteristic fluxes for four-momentum and position. Energy transfer sets
   the label displacement; no endpoint is projected onto a reconstructed chord.
   It matches exact radial and rotating controls, and preserves closure in a
   doublet prefix, but a coarse/fine 90 ps comparison differed by approximately
   **192 confinement reference lengths** in endpoint position. A fine run also
   exposed a nonlinear solver probing an invalid label; a bounded `tanh`
   parameterization allowed its neutral-endpoint continuation to finish.
   Passing the closure gate was insufficient to certify trajectory accuracy.

2. Normalized Hermite ray interpolation with analytic boundary derivatives,
   quadrature of the same stored curve and causal propagation of derivative
   fronts. Manufactured-curve, one-sided jump and checkpoint tests were used.
   Some refinements improved error, but the coarse/fine 90 ps trajectory
   difference was still approximately **7.96 reference lengths**, and sampled
   closure/completion gates failed at tested resolutions. Runs without causal
   derivative-front treatment were substantially worse. This remains an
   experimental method, not a repaired production solver.

These pilots also exposed a diagnostic-guard false positive when the only
remaining interval to a requested observation was a few coordinate ULPs.
An experimental observation-roundoff adjustment was tested, but is not part
of the narrow main-branch patch in this result.

The experimental code and failed-run ledger are kept on
`research/boundary-method-pilot-20261111`. They are not enabled in the main
research solver. Source evolved between development trials, so the ledger
records that they are exploratory comparisons rather than results of a single
frozen revision. No failed case has been relabeled PASS.

## Release decision

| Claim | Result | Evidence / next requirement |
|---|---|---|
| Event traversal and value lookup share the same side convention | PASS | Narrow regression and 26 selected tests |
| Conservative balance laws alone validate an encounter trajectory | FAIL | Coarse/fine discrepancy despite tiny closure |
| Higher-order history pilot is ready for physics production | FAIL | Remaining trajectory and completion failures |
| Original doublet reaches 1000 ps accurately after the narrow fix | UNKNOWN / STOP | No completed long-time counterexample rerun or refined agreement |
| Physical annihilation suppression | UNKNOWN / STOP | Numerical and radiation/quantum matching requirements remain |

The next decisive work is a validated high-accuracy boundary/history scheme
or a targeted late-state reproduction with complete histories. More wall time
on the unchanged failing integrator or loosened closure tolerances would not
establish that. No additional AF job is requested by this change.
