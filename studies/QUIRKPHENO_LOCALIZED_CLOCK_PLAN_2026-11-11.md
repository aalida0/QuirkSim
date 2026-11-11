# Localized proper-clock integration pilot — 2026-11-11

Status: native C++ and synthetic tests pass; full Geant4 build and execution
of the clock update remain pending. This is still a fixed-threshold pion pilot.
Previous evidence of a resolved directional recoil remains recorded separately.

## Implementation

- Replay trial endpoint states without consuming incoming or pushing outgoing
  string vectors. Integrate proper time as three-node Gauss quadrature of m/E
  over lab time. Locate quadrature nodes in the step's curve-length parameter.
- Bracket the 50 ps proper threshold on that computed clock and cap the step
  before the field propagator runs. Update remains the sole operation that
  consumes the accepted incoming history and pushes the outgoing history.
- Measure accepted proper time before Update invalidates history iterators.
  Use that increment in AlongStepDoIt instead of start-of-step m/E times dt.
- Forced sub-navigation steps respect the capped curve length; they may not
  silently take the full available segment. Require moveTinySteps=true.
- Preserve legacy clock behavior as explicitly selected same-input controls.
- Keep the high-precision Geant4 UI conversion enabled, preserving parent mass
  225.1571702763212 GeV and the requested initial momenta.

Root tolerance: 1e-9 ps on the computed clock. Lab quadrature-node inversion
uses 1e-12 ps. The process accepts the numerical root tolerance, and the audit
checks localized event residual magnitude <=2e-9 ps. These are numerical
conditions on the replayed solver path, not rigorous continuum proper-time
error bounds. Field splitting and merged-string approximations still require
trajectory refinement. No event time is snapped to manufacture a pass.

## Local verification

Native C++17 compile with -Wall -Wextra -Werror -pedantic:
- constant-gamma nonlinear curve parameterization and analytic event root;
- relativistic constant-force proper-time integral, including turnaround;
- root locations at multiple starts and lengths;
- invalid brackets, nonmonotone/nonfinite inputs and zero-length clock handling;
- existing pion on-shell/four-momentum and ownership-token tests;
- eight independent Python kinematic fixtures still agree within 1.71e-13 GeV.

Twelve exact upstream preimages checked; new patch checks and applies.
Synthetic ten-case ROOT/log audit passes and rejects corrupted handoff metadata.
Previous uploaded weak archive remains readable and passes its original local
handoff checks. Its rounded benchmark mass remains flagged as mismatched.
No new physical survival, annihilation or continuum causality result is claimed.

## AF suite

Default runner now prepares ten one-event cases at 2 T, Lambda=774 eV, 130 ps:
- stable, endpoint-0 pion +/-y, endpoint-1 pion +y, each at levels 7/8;
- level-8 stable and endpoint-0 +y pion with the legacy clock and identical
  high-precision inputs, isolating clock changes from the input-rounding fix.

Audit measures stable-trajectory changes under the new sampling/clock path,
weak event-time shifts, mass precision, local handoff checks and directional
recoil refinement. The geometric stable control should remain invariant under
clock-only changes; returned data must establish this rather than assumption.

Per-case 120 s plus 5 s kill grace; 600 s build cap; one CPU; no Condor,
set -e or set -u. Stop on first failed case and retain diagnostics.

```bash
cd ~/QuirkSim
git pull --ff-only origin research/quirkpheno-weak-recoil-pilot
batch=/data/aalida/QuirkSim/quirkpheno_clock_$(date +%Y%m%dT%H%M%S)
python3 scripts/af_quirkpheno_weak_pilot.py --output-dir "$batch"
tar --exclude='*/source' --exclude='*/build' -czf ~/quirkpheno_clock.tar.gz -C /data/aalida/QuirkSim "$(basename "$batch")"
```

Upload quirkpheno_clock.tar.gz, including logs on failure. Sampled lifetimes,
leptonic channels and two weak decays remain later tasks. Doublet remains paused.
