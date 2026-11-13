# Physical weak-decay pilot

The runner samples the existing LO Dirac-triplet benchmark: independent exponential
proper thresholds, channels selected by partial widths, and unpolarized rest products.
Neither channel nor lifetime is forced. Free-particle widths are an assumption;
spin correlations, bound-state corrections and higher-order weak corrections remain
outside this calculation. This is fixed-production classical transport, not an LHC
production population or a physical capture/annihilation prediction.

## Default campaign

Neutral mass 225 GeV, Lambda 4642 eV, endpoint beta .4 in the production CM,
common longitudinal boost .25, both endpoints initially charged. Vacuum, uniform
0/2 T, laboratory time limit 1000 ps. Eight independent seeds per field.

Each seed shares exactly the same sampled clocks, channels and rest-product plan
across recoil and forward/backward collinear controls. A control rotates all decay
products together; it conserves local four-momentum but selects an angle rather
than sampling the physical angular PDF. Fields and controls reuse latent events
and are not additional independent statistical events.

Level 8 runs for all seeds; seed 1 also has level 7 counterparts. Four stable
references give 58 native cases total. The two seed-1 zero-field recoil cases run
first. Only a passing long-window conservation preflight releases the remaining
56 cases. The build is shared and direct workers run one CPU per case; no Condor.

The same conservation budgets as the 150 ps Noether pilot remain in force. They
are not increased for the longer duration. Accumulated drift can therefore stop
the preflight. Failures, missing cases and legitimate finite-window censoring all
remain in cohort accounting. Eight events are a pilot, not a precision branching
or tail estimate; rare leptonic modes need not occur and are not injected.

## Encounters and audit

Every endpoint step enters an online asynchronous merge onto common production-CM
time intervals. Minima of piecewise-linear relative position are solved analytically.
Initial and before/after-decay points are included. Duplicate-time parent/daughter
points preserve the momentum/identity discontinuity. Nearby duplicate minima within
1e-9 ps are merged. At most 200000 points per endpoint are buffered; exceeding this
limit aborts explicitly rather than downsampling.

The online calculation is checked against independent offline Python geometry and
known crossings, decay discontinuities and injected buffer/time failures. It records
first and later bracketed minima after both endpoints neutralize, endpoint phase
space and classical orbital vectors. Those vectors are not quantum partial-wave
occupations. Piecewise-linear geometry has no claimed rigorous continuum error
bound; the seed-1 numerical pair reports first-return shifts.

Zero/one-decay histories at the endpoint stop are legitimate censored cases. A
crossed proper threshold without a recorded decay is a failure. A first-return
observable requires both endpoints neutral. Matching is the first minimum after
neutralization, which need not be the same oscillation index across controls.

Direct queue-contour energy/angular checks include emitted-particle birth flux
and the uniform-field axial torque. Full four-momentum and all six orbital
generators are gated at B=0; energy and corrected Jz at Bz. These are classical
sampled conservation checks, not O(hbar) quantum validation. No radiation or capture
process is enabled.

## AF command

On a compute node with four available CPUs:

```bash
cd ~/QuirkSim
git pull --ff-only origin research/quirkpheno-weak-recoil-pilot
out=/data/aalida/QuirkSim/physical_weak_$(date -u +%Y%m%dT%H%M%SZ)
nohup python3 -u scripts/af_physical_weak_ensemble.py \
  --output-dir "$out" --workers 4 >"$out.driver.log" 2>&1 &
echo "$out"
```

Upload <out>_summary.tar.gz. ROOT step trees are disabled. Full logs, event ROOT
files and raw Noether snapshots remain on AF. The compact archive includes the
manifest, audit, input macros/cards/plans, short log tails and online encounter
records, omitting the large per-step conservation streams and source/build tree.

Use --prepare-only to inspect inputs. Use --audit-only with the same output directory
to recover packing without rerunning transport. A preflight failure also produces
a partial summary. No unrun case is promoted to a pass.

The report leaves survival probability null and annihilation uncomputed. Once
paired geometry qualifies, radiation/cooling and capture require explicit physical
models and sensitivity checks.
