# Original-grid targeted AF pilot

Run from the QuirkSim checkout with its existing Python environment active:

```bash
python scripts/af_original_grid_targets.py \
  --output-dir /data/aalida/QuirkSim/grid_targets_$(date +%Y%m%dT%H%M%S) \
  --workers 4
```

The runner builds the pinned public QuirkPheno adapter once, then runs 102
cases with four independent one-CPU processes. It never submits Condor jobs.
Each trajectory has a 600-second default timeout. Completed ROOT files remain
on AF; the streaming audit and packager produce a compact post-neutral archive.
Retain the output directory if a build, transport or audit check fails.

The points are (100 GeV, 2154 eV), (225 GeV, 4642 eV), and (400 GeV, 2154 eV).
Each includes 0/2 T, paired isotropic and forward/backward collinear directions,
four seeds, and selected level-7/8 refinement. Every macro and audit uses its
actual parent mass, neutral mass and string tension.

The C++ adapter is pion-pair conditioned. Proper-time thresholds use the full
LO triplet weak width; the pion-pair fraction and omitted-channel fraction are
recorded for each point. This is a conditional geometry pilot, not a physical
annihilation calculation or an inclusive channel sample.

To check setup without compiling or running trajectories:

```bash
python scripts/af_original_grid_targets.py \
  --output-dir /tmp/quirksim_grid_prepare_check --prepare-only
PYTHONPATH=. python -m pytest \
  research/quirkpheno_recoil/test_af_original_grid_targets.py -q
```

Use a fresh output path for every invocation. Preparation is locally tested;
the new points still require their Geant4 build and transport audit on AF.
