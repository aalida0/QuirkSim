# QuirkSim — standalone charged-to-neutral quirk phenomenology

## Research branch: variational string with internal degrees of freedom

**Branch:** `research/variational-worldsheet` (experimental; do **not**
merge into main or the stable development branch yet).

An independent spatial midpoint discretization of the classical
**Nambu–Goto action** now evolves internal string vertices and their
canonical energy, momentum and angular momentum. It includes a small
explicit unphysical gauge-mass regulator, a uniform solenoidal field,
and a conservative canonical momentum treatment of pion/electron/muon
endpoint decays. It is **not** an independently validated continuum
worldsheet and is not yet able to support paper-quality annihilation
or radiation rates.

A separate encounter model brackets unknown glueball-emission and
annihilation assumptions without claiming physical rates or feeding
radiated energy back into the dynamics. It explicitly conditions on
survival of encounters before a weak decay.

**CRITICAL PHYSICS BLOCKER:** the 260 ps long-time test failed the
spatial mesh-quality gate. A rapidity-coordinate integrator removes the
superluminal numerical integration error, but string segments collapse
below the spatial regularizer and their length ratio exceeds 20 million.
Short-time tests and a 5 ps weak-decay event pass; **multi-oscillation
recoil/annihilation predictions are NOT validated**. See
[GitHub issue #3](https://github.com/aalida0/QuirkSim/issues/3)
for the conservative remeshing/worldsheet convergence acceptance criteria.
Do not use these simulations as paper figures.


**Higher-order worldsheet action validation:** [Nambu–Goto Gauss quadrature,
13/17/21-vertex pre-contact geometric convergence PASS for one triplet
benchmark](studies/NG_GAUSS_QUADRATURE_AUDIT_2026-11-10.md). Full physical
string crossing, weak impulses in this gauge, radiation and annihilation
remain **BLOCKED**. A numerical gate is not a collider prediction.

**Independent conformal-gauge prototype:** moving massive-boundary flux,
causal characteristic histories, exact rotating-string comparison and
four signed radial crossings now have a separate
[derivation](studies/CONFORMAL_BOUNDARY_DERIVATION.md) and
[scoped validation result](studies/CONFORMAL_BOUNDARY_RESULT_2026-11-10.md).
The radial crossing rule is explicitly prescribed in a reduced adapter;
general chart collapse still stops. Electromagnetic boundaries, near-contact
cross-validation and weak recoil in this solver remain pending.

**Mesh-independent first-return validation:** [common total string arclength
convergence and remaining core-scale physics blocker](studies/NG_TOTAL_ARCLENGTH_CONVERGENCE_2026-11-10.md).

**Scientific validation report:**
[variational worldsheet and uncertainty envelope](studies/VARIATIONAL_WORLD_SHEET_PHYSICS_GATES_2026-11-10.md).

Reproduce the action, recoil and transverse-wave checks:

```bash
git switch research/variational-worldsheet
python -m pip install -e '.[test]'
python -m pytest -q
python scripts/variational_gate_report.py
python scripts/validate_transverse_wave.py
python -m quirksim.worldsheet_experiment \
  --seed 25 --duration-ps 5 \
  --output results/variational_weak_seed25.json
```

The clocked example is a **rare, deliberately chosen early-decay
pilot on a pre-stretched string**; neither its contact history nor
its five-picosecond runtime represents an inclusive collider sample.


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

### Physical source of nonradial motion: a uniform solenoid

**New model:** `quirksim/field_dynamics.py` implements the published
Li et al. (2020) transverse-velocity-dependent Nambu-Goto *endpoint-force
approximation* alongside the legacy central-force proxy. A uniform B field
causes (q (\mathbf v \times \mathbf B)) deflection **until each endpoint
weakly decays to its electrically neutral state**.

The benchmark pair starts essentially at a common pp-like production point
with zero initial internal orbital angular momentum; it can carry a
longitudinal pair CM boost (ISR-like kinematic slice). No unphysical
micrometre initial offset is required. This is *not* a full moving
Nambu-Goto worldsheet calculation.

**Three reference commands:**

```bash
# Inspect all field validation tests
python -m pytest -q tests/test_field_dynamics.py

# Paired sample, physical 2T solenoid, published endpoint-force approximation
python scripts/triplet_field_pilot.py \
  --mass 225 --lambda-ev 774 --b-tesla 2 \
  --beta-cm-x 0.25 --force-model ng_endpoint \
  --nevents 20 --output-dir results/field_2T_pilot

# Check the numerical integration and short-distance regulator
python scripts/validate_field_solver.py
```

**Results and caveats:** [2026-11-10 solenoid study](studies/TRIPLET_SOLENOID_NG_ENDPOINT_2026-11-10.md).
The field-on vs field-off pilot and numerical convergence run were completed
on public GitHub Actions. The recoil-induced geometric change persists,
but physical annihilation rates **cannot** be inferred and worldsheet
dynamics remain an open validation gate.

### Triplet nonradial / pre-existing angular momentum sensitivity

The straight-string proxy can be tested away from the artificially radial
production geometry by setting `impact_cm`, interpreted **only** as a
controlled transverse separation between endpoints at simulation start.
For initial back-to-back longitudinal endpoint momenta `+pz, -pz`,

`|L_initial|/hbar = impact_cm * pz / (hbar*c)`.

To scan matched triplet recoil vs no-recoil histories, without touching
other running jobs:

```bash
python scripts/triplet_impact_scan.py \
  --mass 225 --lambda-ev 774 --duration-ps 1000 \
  --impact-um 0,0.01,0.1,1,3,10 \
  --nevents 20 --seed 10 --progress-every 10 \
  --output-dir results/triplet_impact_pilot
```

Outputs include `summary.json`, `summary.csv` and `paired_events.csv`,
with per-offset median post-decay distances, *paired* median differences,
fraction of pairs where recoil increases the distance, contact-reference
fractions, and the initial orbital angular momentum in units of hbar.

**Crucial physics warning:** A pair of quirks from a pointlike pp vertex
does not obtain relative orbital angular momentum merely from a nonzero
opening angle between its initial momenta; in the CM frame that motion is
radial. Offsets of 0.01–10 micrometres are far larger than a hard-production
vertex and are used here exclusively as sensitivity tests. Physically
motivated nonradial motion requires a validated mechanism (magnetic forces,
material interactions, radiative kicks, or a moving-string model). The
`impact_um` scan **cannot** be called a realistic LHC production model.

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

