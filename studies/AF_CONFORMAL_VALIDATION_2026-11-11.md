# AF continuation of the conformal numerical pilots

These jobs finish selected classical validation. They do not provide a
calibrated radiation model or physical annihilation/survival probabilities.
The 1,000-event paired samples from the earlier development branch are a
different calculation and must not be relabeled as conformal-solver results.

## Submit the interrupted triplet seed 1 pilot

Use the existing AF checkout and standalone Python virtual environment:

```bash
cd ~/QuirkSim
git fetch origin
git switch research/variational-worldsheet
git pull --ff-only origin research/variational-worldsheet
source .venv/bin/activate
python -m pip install -e '.[test]'
mkdir -p results
python -u scripts/validate_conformal_evolution.py \
  --mode recoil --model triplet --seed 1 --duration-ps 1000 \
  --max-step-ps 0.025 --ray-angle 0.0001953125 \
  --causal-fraction 0.000390625 --max-refinements 0 \
  --output results/conformal_triplet_seed1_strict.json
```

If the research branch is not present locally, `git switch` should create it
from the uniquely named remote tracking branch after fetching. The helper
requires this branch and refuses tracked uncommitted changes. No Athena or
Key4hep environment is required. Run this single pilot directly on an AF compute node. Condor is reserved
for the parallel convergence batch below. Seed 1 has since been reported
PASS by the user; it need not be rerun just to recover the interruption.

The parallel submission helper prints its exact pinned commit, manifest and result directory.
It archives committed source under `~/QuirkSim/condor/`, with Condor logs
there, and writes physics output under
`/data/aalida/QuirkSim/conformal_validation/`. Every invocation has a unique
directory. It records the Python executable and dependency versions; the
shared `.venv` must remain available and unchanged while jobs run. Source and manifest are transferred into Condor scratch; output is written
there and transferred back under `~/QuirkSim/condor/RUN/returned/CASE/result/`.
All Condor working directories, transfer destinations and logs are under home;
AF forbids submitting jobs with an initial directory under `/data`. The
`collect` command copies completed results to `/data` after transfer. Worker
code does not read home source paths or write `/data` directly. The Python
runtime still uses the absolute shared AF `.venv` path. The wrapper checks
that it exists; this change does not make that runtime portable to workers
without AF home access. A missing runtime is reported in `bootstrap_error.txt`.

The initial submission used `should_transfer_files = NO`, which imposed a
FileSystemDomain match and rejected all advertised slots (cluster 4252968).
The corrected helper uses `YES`, explicit inputs and output directory
transfer. Do not merely remove the scheduler requirement from the old job.

Each job requests **one CPU, 4 GB memory and 2 GB scratch disk**; BLAS thread
counts are fixed to one. A job is removed after eight running hours. Site
policy can impose a shorter limit. An interrupted or removed job remains
INCOMPLETE, never PASS. Logs flush at each solver progress report.

`remaining` runs triplet recoil, seed 1, for 1,000 ps with step 0.025 ps,
ray angle 0.0001953125 and causal fraction 0.000390625. Automatic refinement
is disabled so a failed resolution is preserved explicitly.

## Subsequent convergence batch

After checking the remaining pilot's result:

```bash
python scripts/af_conformal_validation.py prepare --suite convergence --submit
```

This submits eight independent jobs: doublet/triplet, seeds 1/25 and two
resolutions. The first resolution is the strict configuration above; the
second halves all three integration controls. The fixed benchmarks are
225 GeV/774 eV for the triplet and 150 GeV/2154 eV for the doublet. Both use
the same stretched initial string and boost/field settings as the existing
conformal pilot. This is not a production-kinematics ensemble.

The overlapping seed 1 strict run intentionally provides a fresh result
within the same pinned-source batch. There is no automatic chain that
submits the larger batch after a failed first run.

## Monitor and collect

```bash
condor_q
```

The submission output gives exact paths; inspect `scheduler.log` and each
job's `.out`/`.err` under its home directory. Under the printed result
directory each case contains a transferred `result/` folder with:

- `provenance.json`: source commit, integration settings and command;
- `progress.log`: flushed progress and solver errors;
- `evolution.json`: complete result, including a numerical FAIL when present;
- `completion.json`: process return code; absent for interruption/removal.

Collect and summarize using the **exact manifest path printed at preparation**:

```bash
python scripts/af_conformal_validation.py collect ~/QuirkSim/condor/EXACT_RUN_DIRECTORY/manifest.json
```

This saves `summary.json` and exits nonzero if any case failed or is
incomplete. Return the summary plus the `evolution.json` files for comparing
resolutions. The gate checks closure below 0.1 confinement-reference lengths,
energy residual below 1e-7 GeV, momentum residual below 5e-6 GeV and an accepted
termination. Two independent PASS flags alone do not prove convergence:
matching times, endpoint trajectories and weak events must also agree at the
scale used by the intended observable. No physical survival flag is enabled.

## Helper validation

Before publication of this helper, a temporary git fixture verified source
archiving, manifest generation, worker execution with the actual solver for
0.1 ps, result collection and summary. The smoke calculation passed its
numerical gate. A nonzero worker return code could not become a summary
PASS; rerunning over existing result files was refused. The eight-case matrix
was checked. These are submission-helper checks, not completion of the AF
1,000 ps pilots. Actual scheduling and compute-node runtime require AF.

The worker uses explicit error checks and no `set -e` or `set -u`.
Collection can be repeated; incomplete cases remain INCOMPLETE and differing
previously collected results are never silently overwritten.
