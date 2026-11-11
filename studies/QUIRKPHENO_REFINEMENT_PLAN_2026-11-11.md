# Stable triplet-scale QuirkPheno refinement suite — 2026-11-11

Status: code prepared and locally checked; new AF simulations not run in this session.
Upstream remains pinned to d91cb27ed3786bc63900faad892513841b5b6ed4.
Mass 225 GeV is a stable transport control, not the loop-split charged parent.
Lambda=774 eV; beta endpoints=.4 in pair CM; beta CM z=.25.
Doublet and weak-recoil production remain paused.

## Cases

Eight sequential single-event runs, 0/2 T and levels 1–4:

| Level | maxBoost | maxMerge (mm) |
|---|---|---|
| 1 | 0.005 | 5e-9 |
| 2 | 0.0025 | 2.5e-9 |
| 3 | 0.00125 | 1.25e-9 |
| 4 | 0.000625 | 6.25e-10 |

Nominal end 130 ps. In zero field the analytic first turnaround is about
111.43 ps and first re-encounter about 222.86 ps. Thus this test targets
turnaround, not crossing/quantum matching. The magnetic case must be checked
from its own output, not assigned those exact analytic times.
Vacuum, every step saved, crossing restart disabled, one CPU.
Both controls are refined together: this tests combined sensitivity and
does not identify a unique error source or convergence order.
Same initial-condition mismatch with the stretched conformal pilot remains.

Per-case timeout 120 seconds plus 5-second kill grace. Eight cases cap at
roughly 17 minutes plus build (capped at 10 minutes) and setup. These are
limits, not runtime predictions. No Condor or set -e/-u. Existing 20 ps
smoke suite stays available as the default.

## Audit changes and local evidence

The audit accepts available levels, compares adjacent resolutions through
the common requested/covered interval without extrapolation, and records first
saved px-sign-change brackets. It does not infer a precise root from sparse
saved observations or equate a sign change with quantum matching.
The zero-field radial analytic check extends through turnaround but rejects
data at/after the first analytic re-encounter. No magnetic interpolation
error certificate or total string conservation claim is supplied.

Smoke gates now use actual execution status, macro hashes, one event,
two endpoints, finite/increasing records, row counts, terminal kTime reason
and requested duration coverage. Missing ROOT files are explicit failures.
A smoke pass does not mean accuracy passed. GNU timeout exit codes 124/137
are marked as timeouts; an incomplete runner exits nonzero and retains logs.

Local checks completed:
- eight generated cards/macros; correct 130 ps conversion and levels;
- bytecode compilation for runner and audit;
- previous uploaded archive: all four numerical resolution comparisons
  reproduce to 1e-18 mm;
- timeout exit 124 recognized;
- analytic turnaround/re-encounter ordering for the selected duration.

No C++ changes were made. Geant4 runtime remains unavailable locally.

## AF command, directly on a compute node

```bash
cd ~/QuirkSim
git pull --ff-only origin research/variational-worldsheet

batch=/data/aalida/QuirkSim/quirkpheno_refinement_$(date +%Y%m%dT%H%M%S)
python3 scripts/af_quirkpheno_transport.py --suite refinement --output-dir "$batch"

tar -czf ~/quirkpheno_refinement.tar.gz -C /data/aalida/QuirkSim "$(basename "$batch")"
```

Upload the archive including logs, manifest and ROOT outputs. Keep any timed-out
or failed cases. No recoil adapter should be promoted before reviewing the
magnetic refinement and turnaround evidence.
