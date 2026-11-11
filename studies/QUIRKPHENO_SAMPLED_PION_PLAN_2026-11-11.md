# Small sampled-lifetime pion-channel pilot — 2026-11-11

Status: prepared; native and synthetic checks pass; new Geant4 execution pending.
The preceding ten-case localized-clock AF pilot passed all local checks and
selected trajectory comparison gates. No new simulation was launched here.

## Scope and sampling

Fixed initially charged conjugate pair, neutral mass 225 GeV, charged mass
225.1571702763212 GeV, Lambda=774 eV, Bz=2 T, beta=.4 in pair CM,
CM beta_z=.25, point production and nominal observation through 130 ps.

Each endpoint receives an independent exponential proper threshold with the
existing nominal LO total weak lifetime 258.5200888747253 ps. The rest-frame
pion direction is isotropic. Channels are forced to pion, not sampled from
all partial widths. The reference pion branching 0.976048 is recorded but
must not be mistaken for an applied inclusive normalization. No physical
survival, production rate, angular distribution estimate or efficiency is
certified by eight seeds. Free-particle weak rates omit bound-state/IC and
higher-order corrections. Reference and backend pion masses differ by 0.29 keV;
the earlier result report records this convention.

Predeclared complete seeds 1-8; no cherry picking or filtering by decay count.
Each has a physical antipodal-direction companion with identical proper-time
thresholds. Subsequent lab event times may differ physically after the first
recoil; do not require the second lab time to match. Both endpoint slots keep
sign-preserving ownership across charged-to-neutral transitions.

Main eight pairs run at level 8. Seed 1's pair also runs at level 7. Two stable
references bring the total to 20 one-event cases. The seed block includes
thresholds capable of zero, one and two decays in this finite window; actual
counts must come from returned trajectories. Lifetimes beyond the observation
window are censored, not failures or annihilation.

## Implementation and evidence

Versioned two-row plans contain proper thresholds (ps) and unit rest directions.
They are loaded immutably once per worker in a one-event process. A shared
endpoint-specific threshold function is used by transport root localization
and the weak process; both charged processes can be enabled. Neutral species
have no weak process. Parent/daughter and partner/string preservation checks
from the deterministic pilot remain active. No string restart is introduced.

Local checks:
- strict C++17 native core, clock and plan parser tests;
- invalid/misordered plans, lifetimes, directions and trailing data rejected;
- exact Python-to-C++ roundtrip for all 18 sampled plan files;
- complete seed block, shared thresholds and antipodal unit directions;
- synthetic 20-case ROOT audit exercises 0/1/2 decays and corruption rejection;
- prior deterministic ten-case synthetic audit still passes;
- patch checks/applies against the same 12 exact upstream preimages.

Synthetic fixtures are not Geant4 validation. Seed-1 paired trace comparison
and causal checks remain required after execution. Leptonic channels, production
kinematics, matching, radiation and annihilation remain beyond this pilot.
Doublet stays paused.

## Direct AF compute-node run

Each case has a 120 s cap plus 5 s kill grace; build cap 600 s. Stop on first
failed case; retain logs. No Condor or shell set -e/-u.

```bash
cd ~/QuirkSim
git pull --ff-only origin research/quirkpheno-weak-recoil-pilot
batch=/data/aalida/QuirkSim/quirkpheno_sampled_pion_$(date +%Y%m%dT%H%M%S)
python3 scripts/af_quirkpheno_weak_pilot.py --sampled-pion-seeds 8 --output-dir "$batch"
tar --exclude='*/source' --exclude='*/build' -czf ~/quirkpheno_sampled_pion.tar.gz -C /data/aalida/QuirkSim "$(basename "$batch")"
```

Upload the archive with manifest, plans, macros, ROOT files and logs.
