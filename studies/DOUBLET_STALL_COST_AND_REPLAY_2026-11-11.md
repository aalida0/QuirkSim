# Three-hour doublet diagnostic: precision limit and bounded replay

Evidence: user-provided AF stdout for the 330 ps doublet seed 1 diagnostic
at the instrumented `4d983cb` code. The quoted three-hour duration is
user-reported. This session did not receive a new complete evolution JSON or
the internal characteristic histories from that run. The earlier eight-case
archive remains the independently inspected full-output evidence.

## What the diagnostic establishes

The run stops at 315.081015984348 ps. In its **last advance call**, not the
whole job, there were 162655 attempted steps, 2529693 rejected ray trials and
30704 event-root truncations. The final proposed step was
8.282815571546829e-14 ps, whereas both outgoing-coordinate ULPs were
5.684341886080802e-14 ps. The step is only 1.45713 coordinate ULPs. One proposed
outgoing coordinate increments by zero, the other by exactly one ULP.
The last ray trial required 27 angle-based halvings and had no event roots.
Thus the terminal failure is an inability to represent a strictly increasing
outgoing null coordinate after extreme local resolution work. Raising the
200000-step budget will not resolve that floating-point limit.

This is not a physical chart-collapse determination: the energy-coordinate
interval is 0.02845954 ps versus the configured 1e-8 ps cutoff. Total energy
residual remains tiny (3.43e-11 GeV), but maximum geometric closure is 1.14815
reference lengths, exceeding the 0.1 threshold. The run was already unsuitable
for core-scale inference before it reached the terminal precision failure.

The output does **not** distinguish a cascade of physically blueshifted ray
features from numerical interpolation/reflection defects. The final incoming
jump cap is not an event root of that last trial. It would be unjustified to
attribute the entire stall to a repeated zero-distance event root on these
numbers alone. Neither annihilation nor loss of physical survival is inferred.

The earlier instruction to run the full 330 ps prefix lacked an adequate cost
limit. Do not repeat that unbounded diagnostic.

## New opt-in diagnostic controls

`validate_conformal_evolution.py --bounded-diagnostic` provides:

- 120 seconds of integration time, shared across its observation intervals;
- at most 50000 rejected ray trials per `advance` call;
- a conservative stop when a proposed trial step is no larger than 16 ULPs
  at the current time/outgoing coordinates;
- a compressed `.npz` checkpoint next to the JSON result.

These limits return explicit numerical/resource STOP reasons. They do not
accept rejected trials, relax the closure threshold, move endpoints, change
the force law or create a quantum matching event. They are opt-in; default
integration controls remain unchanged. The 120-second bound covers integration,
not imports, compression and final serialization; a single expensive operation
can cross the deadline before the next check. The 16-ULP threshold is an early
diagnostic precaution, not a proof of sufficient floating-point accuracy.

The checkpoint preserves endpoint state, masses/charges, time/force/field,
both ray coordinate and direction histories, their already accumulated
integrals, discontinuity lists and weak/reflection records. Loading uses
`allow_pickle=False` and does not recompute histories or sums. This is a solver
checkpoint, not a complete event-scheduler checkpoint. The standalone replay
therefore requires both endpoints already neutral, so no future weak decay
is silently omitted. Checkpoints generated earlier in the weak evolution
cannot be replayed with that standalone script.

`scripts/replay_conformal_checkpoint.py` replays a small post-decay interval
under a 30-second integration budget and a 5000-step budget. It is a classical
kernel diagnostic only: it does not apply a quantum matching surface or claim
physical survival. Every replay should start afresh from the same checkpoint.

The previously completed three-hour JSON does not contain the full histories,
so it cannot be turned into this checkpoint after the fact. No further long
AF calculation is requested now. A future bounded run may stop before the
original failure; its checkpoint is for local method comparisons, not evidence
that the 315 ps stall has been reproduced under the new limits.

## Validation and next repair criterion

Tests require an exact checkpoint round trip, including a pair of weak-induced
ray jumps and accumulated primitives, and identical subsequent accepted state.
Resource/resolution guard tests require no unaccepted endpoint mutation and
explicit failure rather than a quantum-survival result. Existing diagnostic
tests remain applicable. The focused set passes 24 tests, including six checkpoint/diagnostic tests.
The separate replay CLI smoke test reaches its short target from a neutral
checkpoint. The full project suite and the failed 315 ps case were not rerun.
A bounded CLI smoke run stopped on its 120-second integration limit at
39.0745 ps and wrote a checkpoint. Only one weak decay had occurred, so
this particular checkpoint is intentionally rejected by neutral-only replay.
It verifies the resource guard, not the late stall. No claim of a repaired
doublet is made. Results are recorded in
`data/conformal_checkpoint_limits_verification_2026-11-11.json`.

The next numerical repair must explain and control the fine ray structure,
preferably through a comparison against the segment-wise propagator identified
in the literature audit. It must recover matched-time closure and trajectory
convergence without projecting positions, merging away unresolved physical
features, increasing acceptance thresholds, or reclassifying failures as contact.
Checkpoint and budget infrastructure is not itself that repair.
