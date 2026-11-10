# QuirkSim — standalone charged-to-neutral quirk phenomenology

Independent, publicly reproducible project. **No ATLAS-internal code.**

## Current prototype

The simulator includes Dirac vectorlike SU(2)L doublet (Y=-1/2) and triplet (Y=0)
mass splittings (one-loop approximation), leading pion/electron/muon branching
fractions, independent proper-time clocks for two charged endpoints, exact
event-by-event four-momentum conservation in two- and three-body weak decays,
and numerical evolution with a **central straight-string potential proxy**.

The 3-body rest-frame sampler implements the leading unpolarized
vector-fermion-current x V-A leptonic-current angular and invariant-mass
distribution, with widths separately fixed using published LO formulas.

**This is NOT a full Nambu–Goto worldsheet solver.** The string is represented
by a regulated central potential; no internal worldsheet dynamics, gluonic
radiation, pion-medium effects, annihilation, magnetic fields, detector matter,
or realistic pp production kinematics are yet implemented. No conclusions
about physical annihilation probabilities or LHC reach are supported.

## Install

```bash
git clone https://github.com/aalida0/QuirkSim.git
cd QuirkSim
git switch development/standalone-simulator
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest -q
```

You do **not** need Key4hep, Athena or Geant4. If using CVMFS/Key4hep anyway,
prefer `python -m pytest` to avoid a mismatched external `pytest` executable.

## Single-event examples

```bash
python -m quirksim.cli --model triplet --mass 225 --lambda-ev 774 \
    --mode recoil --duration-ps 1000 --seed 42 --trace \
    --output results/triplet_recoil.json

python -m quirksim.cli --model doublet --mass 150 --lambda-ev 2154 \
    --mode recoil --duration-ps 1000 --seed 1 \
    --output results/doublet_leptonic.json
```

The second example explicitly exercises an electron-neutrino decay, which is
now fully evolved rather than terminating as `leptonic_unresolved`.

## Matched paired ensembles

```bash
python scripts/paired_compare.py --model triplet --mass 225 --lambda-ev 774 \
    --duration-ps 1000 --seed 10 --nevents 100 \
    --output-dir results/paired_triplet

python scripts/paired_compare.py --model doublet --mass 150 --lambda-ev 2154 \
    --duration-ps 1000 --seed 10 --nevents 100 \
    --output-dir results/paired_doublet
```

### Performance and interrupted-run safety

The three-way ensemble uses three independent adaptive SciPy ODE integrations
per seed, and large samples can take many minutes. The script now prints
progress (completed seeds, elapsed time, average seconds/seed, estimated
remaining minutes) every 10 seeds, and writes atomic partial `events.csv`
and `pairs.csv` snapshots every 25 completed seeds. `summary.json` is
written only after the complete run. Partial snapshots are **not automatically
resumed**; use distinct seed ranges when splitting a large run.

To omit the redundant stable reference and run just the matched recoil
and no-recoil cases, pass `--skip-stable` (two ODE simulations per seed
rather than three). This reduces work by roughly one third; it is not
a faster integrator.

On a shared login node, run only small pilots. For example:

```bash
python scripts/paired_compare.py --model triplet --mass 225 --lambda-ev 774 \
    --duration-ps 500 --nevents 20 --skip-stable \
    --progress-every 5 --checkpoint-every 5 \
    --output-dir results/pilot20
```

Submit larger runs to a compute/batch node (e.g. HTCondor); do not launch
large multi-core jobs on login nodes. Physics comparability requires keeping
the same observation window when comparing models.

Each run outputs `events.csv` (three modes), `pairs.csv` (matched decay
histories) and `summary.json`. Each endpoint reuses its sampled *proper*
lifetime, chosen decay channel, and rest-frame random seed in both modes.
**Lab-frame decay times need not match** after the first recoil.

- `stable`: no weak decays, charged endpoints throughout.
- `no_recoil`: both endpoints weak-decay but momentum is artificially unchanged
  (unphysical control; does not conserve whole-system energy at transitions).
- `recoil`: explicit two-/three-body decay products, momentum-conserving
  endpoint updates and continuation through both decays.

### Observables, limitations and interpretation

`minimum_separation_including_production_cm` contains the trivial near-contact
production geometry and **must not** be interpreted as post-decay contact.
`turning_points` contain local minima of endpoint separation, found by event
root-finding for d|r|²/dt=0. `closest_post_decay_turning_point_cm` is measured
only after the first weak decay. Values are `null` if no such minimum occurs
within the observation window (not zero).

The geometric core radius `hbar*c/Lambda_IC` is a **reference length**, not an
annihilation criterion. Encounters closer than this radius are *geometric
contacts*, not physical infracolor-emission or quirkonium-annihilation rates.
Ensemble fractions are conditional on producing an *initially charged*
pair with fixed input beta, not inclusive pp production predictions.

**Known important bias:** radial production geometry guarantees perfect
zero-impact-parameter returns in stable/no-recoil controls. A nonzero impact
parameter or realistic 3-D momentum distribution and moving-string
worldsheet are mandatory before asserting a genuine recoil effect.
The no-recoil control is deliberately unphysical, not a competing model.

## Required validation before paper-quality phenomenology

1. Independent physics audit of pion and leptonic normalization and
   electroweak mass-splitting input schemes, with quoted uncertainties.
2. Systematic solver convergence for nonzero impact parameter and regulator.
3. Correct moving Nambu–Goto worldsheet dynamics, including momentum transfer
   into/out of string degrees of freedom.
4. Radiation and annihilation models with explicit uncertainty variations,
   not identification with a geometric encounter.
5. Realistic parton-level pp production via an independent UFO and MG5.

## Core literature

- Thomas & Wells, *Massive Doublet Leptons* (1998), hep-ph/9804359.
- Kang & Luty, *Macroscopic Strings and Quirks at Colliders* (2009), 0805.4642.
- Martin, *Quirks in Supersymmetry with Gauge Coupling Unification* (2011), 1012.2072.
- *The Quirk Trajectory* (2020), 2002.07503.

Development should remain on the feature branch until independently validated.
