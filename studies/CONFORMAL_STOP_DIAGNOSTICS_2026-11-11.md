# Targeted doublet diagnostics

Instrumentation exposes the last failed `advance` call in the evolution JSON:
reason, current and target times, attempted steps, per-call budget, energy
coordinate interval, ray-history sizes, last/null proposed coordinates and
their increments/ULPs, last trial step, angle/rejection count and event caps
and roots when available. Aggregate rejected ray trials, root truncations and
minimum accepted step distinguish budget cost from the last single trial.
At a budget stop, the trial fields describe the preceding accepted step;
the explicit `trial_context` label prevents treating its zero outgoing-label
increments as a failed proposal.

The additions do not change the integrator equations or default stepping.
They do not repair the doublet or enable a physical survival prediction.
`--max-steps-per-advance` exposes the existing budget (default remains 200000).
An optional `--stop-on-closure-failure` stops at the first saved observation
violating the unchanged 0.1-reference gate. Its reason is numerical
`STOP_SAMPLED_CLOSURE_TOLERANCE`, never quantum matching. With this flag off,
the original failed trajectory can continue to its floating-point/budget
stop for diagnosis. Observation checkpoints are not dense output.

## Local reproduction and tests

The doublet seed 1, level 0 prefix reproduced a first sampled closure failure
at 78.23578739115479 ps, the same observation time as the uploaded AF archive.
Closure was 1.3272624656713997e-9 cm = 0.14488254 reference lengths; energy
residual was -3.979e-13 GeV. This is deliberately a FAIL result, saved in the
verification JSON. It identifies early inaccuracy, not the later stall's
root cause. Python/platform differences need not give bit-identical roundoff
against the AF archive.

Validation includes the existing characteristics, impulse and matching tests
and three new diagnostic tests (21 focused tests pass). A budget-limited
trajectory is checked against the corresponding successfully advanced state;
diagnostic serialization and clearing on a subsequent successful call are
checked. A numerical stop never supplies a physical survival probability.
An original-versus-instrumented serialized 15 ps triplet seed 25 comparison
with two weak decays checks that default numerical outputs remain identical
after excluding added metadata. No new 1,000 ps validation is claimed and
the full project test suite was not rerun in this work session.

## Next AF diagnostic: one job directly on the compute node

Do not submit another eight-job batch. To capture the late seed 1 doublet
stall, update the research branch, activate the existing `.venv`, then run:

```bash
mkdir -p results
python -u scripts/validate_conformal_evolution.py \
  --model doublet --mode recoil --seed 1 --duration-ps 330 \
  --max-step-ps 0.025 --ray-angle 0.0001953125 \
  --causal-fraction 0.000390625 --max-refinements 0 \
  --output results/doublet_seed1_stall_diagnostic.json
```

This intentionally omits the early closure-stop flag to expose the original
late failure. A nonzero exit is expected if the failure reproduces. It is
not a success run or a result for physics inference. No `set -e`/`set -u`
shell flags or Condor submission are required.

For the faster early-accuracy reproducer instead, use duration 90 and add
`--stop-on-closure-failure`; that check has already been performed locally.

Do not increase budgets, merge near-null samples, project endpoint positions
onto the chord, or alter the closure threshold on the basis of the current
evidence. The [reuse search](QUIRK_REUSE_SEARCH_2026-11-11.md) identifies a
segment-wise hyperbolic propagator and massive-boundary characteristic
literature as candidates for the next independent formulation comparison.
