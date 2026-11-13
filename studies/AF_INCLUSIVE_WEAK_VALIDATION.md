# Inclusive weak-decay adapter: small native AF validation

This is a validation fixture, not a physical decay-time ensemble or an annihilation result. It extends the existing pinned public QuirkPheno handoff to electron/neutrino and muon/neutrino channels while retaining the pion-only interface for previous runs.

## What runs

40 independent one-CPU runs, four concurrent workers by default. No Condor submission, shell `set -e` or `set -u` is used by the driver. Build once; 600-second transport limit per case.

- Neutral triplet mass 225 GeV, Lambda=4642 eV; charged mass from the existing QuirkSim LO benchmark.
- Pair CM endpoint beta=0.4 and common z boost=0.25.
- Uniform 0 and 2 T in vacuum; 500 ps lab-time limit.
- Pion, electron+neutrino, and muon+neutrino channels explicitly forced at the two endpoints. Proper thresholds 50 and 100 ps ensure every channel is exercised without waiting for rare decays.
- Each weak channel has isotropic and physical forward/backward collinear controls, with identical input rest-product invariants. Every mode/field is repeated at numerical levels 7 and 8.
- Four stable-parent runs complete the suite.

There are no collider event weights, inclusive branching fractions, fitted radiation, capture probabilities or detector efficiencies in this fixture. Its clocks do not represent the physical exponential lifetime distribution. The ordinary channel sampler remains available in `quirksim/weak_plan.py` without forced arguments.

## Run on an AF compute node with four CPUs available

Use the branch containing this runner. Keep sources in home and outputs in data. Use the existing Python/LCG-compatible environment; the runner sources the same LCG_110 setup as previous successful QuirkPheno campaigns.

```bash
cd ~/QuirkSim
git fetch origin
git switch research/quirkpheno-weak-recoil-pilot
git pull --ff-only origin research/quirkpheno-weak-recoil-pilot
python3 -m pip install uproot
out=/data/aalida/QuirkSim/inclusive_weak_$(date -u +%Y%m%dT%H%M%SZ)
nohup python3 -u scripts/af_validate_inclusive_weak.py \
  --output-dir "$out" --workers 4 >"${out}.driver.log" 2>&1 &
echo "$out"
```

Monitor with `tail -f "${out}.driver.log"`. The build runs before any transport case. If CVMFS or compilation fails, retain the logs; this is not a physical failure. Do not submit a single serial Condor job.

At completion, return **`<output-directory>_compact.tar.gz`**. ROOT files stay on AF. The archive contains input cards, generic weak plans, vertex/handoff checks, selected neutral return windows and build logs. It does not include the large raw ROOT files.

If the connection drops during auditing/packing, recover without transport:

```bash
python3 scripts/af_validate_inclusive_weak.py --output-dir "$out" --audit-only
```

The output directory must be new for preparation/execution. `--audit-only` is the explicit reuse path. To inspect all cards without Geant4:

```bash
python3 scripts/af_validate_inclusive_weak.py --output-dir /tmp/inclusive_weak_prepare --prepare-only
```

## Acceptance and limitations

1. Compare input hashes, process return codes and ROOT presence. Require both forced decays, their matching proper thresholds and neutral masses.
2. Independently recompute local four-momentum conservation at 1e-6 MeV tolerance and clock residuals at 2e-9 ps. Product charges/flavors and on-shell masses are checked by the native plan reader.
3. Preserve endpoint string histories through the existing fingerprint/handoff mechanism; reject reported orphan endpoints.
4. Compute first post-neutral common-CM return estimates using bounded streaming ROOT windows. Store level-7/8 shifts and explicit failed/censored cases. Piecewise-linear minima are estimates, not formal transport-error bounds or capture labels. Compare the two levels; do not silently turn their agreement into continuum certification.

Standard product masses are explicitly stored in the plan. Geant4's particle table can differ by a fraction of a keV. The inclusive process permits at most 1 keV disagreement, assigns the planned dynamic mass/kinetic energy, and checks the actual resulting momentum before inserting tracks. This avoids an artificial local conservation error caused by mixing mass tables. The legacy pion path is unchanged.

Whole-string Noether conservation, arbitrary transverse production boosts, initial charged-neutral topologies, full mass-grid coverage and quantum annihilation are separate checks. This gun fixture does not certify them.

## Local verification performed

18 focused tests passed: new generic native parsing/boost/rotation checks, preparation/pairing, malformed-vertex rejection, analytic common-time minima, streaming ROOT round trip, existing complete-decay tests and original-grid preparation regression. The old independent C++ clock/handoff/direction checks also pass. The pinned upstream patch applies cleanly at `d91cb27ed3786bc63900faad892513841b5b6ed4`.

The Geant4 compilation and 40 transport runs are **not claimed performed locally**: this environment has no AF/CVMFS installation. The AF run is the next execution gate.
