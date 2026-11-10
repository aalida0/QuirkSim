# Magnetic-field and straight-string endpoint-force validation — 2026-11-10

**Project:** public standalone QuirkSim, completely independent of Jackson's
private ATLAS analysis. **Status: PHYSICS_CONDITIONAL / NO FULL-WORLDSHEET
VALIDATION.** The results below are deterministic reproducible pilots, not
LHC observables or predicted annihilation probabilities.

## Physical mechanism and method

Unlike the previous unphysical initial transverse-offset scan, both charged
quirks are now created at essentially **one space point** (regulator separation
1e-7 cm, aligned with the relative momentum), with exactly zero initial
internal orbital angular momentum. They have opposite electric charges (+1,-1)
and equal masses. The pair has relative CM beta = 0.4 and a common lab-frame
boost beta_CM,x=0.25 along their production axis, which could be obtained in
pp production with a recoiling ISR system. This is a controlled kinematic
slice, not an MG5 event distribution.

Inside a uniform solenoid B=(0,0,Bz) Tesla we include:

dp_i/dt = F_string,i + q_i K (v_i/c x B),
K = 0.299792458 GeV/(T*m) * 0.000299792458 m/ps.

**Charge is set to zero for each endpoint after its weak decay.** We
explicitly simulate pion/electron/muon channels and preserve the sampled
proper-time decay clocks in both physical-recoil and artificial no-recoil
controls.

We compare two distinct string-force choices:

1. `central_proxy`: smooth central force toward the other endpoint,
   used in QuirkSim's original pilot.
2. `ng_endpoint`: *published* approximate straight-string endpoint force:
   F_s = -sigma sqrt(1-v_perp²) s_hat -
         sigma v_parallel v_perp / sqrt(1-v_perp²),
   evaluated with an instantaneous lab-frame separation direction s_hat.
   This follows Li et al., EPJC 80 (2020) 651, Eq. (2.2),
   **but is NOT a complete dynamically evolving Nambu-Goto string**.
   String curvature, its internal momentum and angular momentum, and
   finite propagation of endpoint disturbances are not represented.

For both, the near-coincidence force is smoothly regulated to prevent an
unphysical singular unit vector. The regulator is varied in numerical tests.

### Primary references

- Li, Li, Pei, Zhang, *The quirk trajectory*, EPJC 80, 651 (2020):
  https://link.springer.com/article/10.1140/epjc/s10052-020-8209-y
  Eqs. (2.1)-(2.3), uniform magnetic fields Sect. 2.2 and the caveat on
  the changing CM frame in Sect. 2.
- ATLAS public magnet description: https://www.atlas.cern/Discover/Detector/Magnet-System
  (nominal solenoid field **2 Tesla**).
- Kang & Luty, *Macroscopic Strings and Quirks at Colliders*, arXiv:0805.4642.

## Simulations

Benchmark: Dirac SU(2)L triplet Y=0 with initial charged +/- endpoints,
m_neutral=225 GeV, Lambda_IC=774 eV, beta_rel=0.4, beta_CM,x=0.25,
1000 ps, 50 independently sampled pairs per model, seeds 10-59,
paired no-recoil/recoil trajectories per seed.

Run: https://github.com/aalida0/QuirkSim/actions/runs/42835314370
(Four successful jobs and their CSV/JSON artifacts.)

All distances in micrometres. The 'paired increase' is the **median of
eventwise recoil-minus-no-recoil closest-approach differences**, rather than
the difference of the separate medians.

| String model | Bz (T) | median pre-decay internal \|L\| / hbar | median first recoil \|Delta L\|/hbar | median no-recoil rmin | median recoil rmin | median paired rmin increase | recoil increases rmin |
|---|---:|---:|---:|---:|---:|---:|---:|
| central_proxy | 0 | 0 | 1.419e12 | ~0 | 1.016 | +1.016 | 50/50 |
| central_proxy | 2 | 1.513e11 | 1.419e12 | 0.060 | 0.887 | +0.784 | 48/50 |
| ng_endpoint | 0 | 0 | 1.419e12 | ~0 | 3.046 | +3.046 | 50/50 |
| ng_endpoint | 2 | 1.653e11 | 1.419e12 | 0.401 | 3.068 | +2.774 | 49/50 |

*The 'pre-decay internal L' is the norm of a kinematic relative orbital
angular-momentum proxy (r1-r2) x (p1-p2)/2 at the first weak decay.
It is not the complete field+string angular momentum.*

In the 2T `ng_endpoint` pilot, field-generated L before weak decay is
~12% of the median weak-recoil first-kick angular momentum. Crucially,
**the field removes the artificially exact rmin=0 no-recoil control**.
Recoil still gives a positive paired increase in 49/50 events.

### Geometric core contact: DO NOT interpret as annihilation

With the naive reference radius r_core=hbar*c/Lambda_IC, the
`ng_endpoint` 2T sample gives 3/50 no-recoil events and 0/50
recoil events with at least one post-decay geometric core encounter.
The central-force 2T sample gives 0/50 versus 0/50. Both findings
are small-sample and strongly model-dependent. They **do not** imply
any experimentally meaningful survival/annihilation difference.
Events may undergo string interactions and annihilate on earlier
returns *before* their simulated weak decays; those processes are
currently absent. A simplistic core-distance crossing is not a
physical annihilation model.

## Validation gates and major remaining blockers

**Numerical sanity: PASS for the tested subset.**

- Published EoM endpoint-force transverse and radial projections tested
  algebraically against Eq. (2.2).
- With B=0 and radial initial motion, both string-force implementations
  reproduce the analytical first-return time to <0.1 ps and retain the
  simple endpoint+string energy diagnostic to <1e-4 GeV.
- No B and initial coincident radial pair gives zero internal orbital L
  even with a nonzero pair boost. A 2T field with asymmetric laboratory
  endpoint velocities produces nonzero internal orbital L; reversing B
  reverses its direction (magnitudes agree).
- Correct pion, electron and muon weak decays, endpoint charge switch-off,
  paired proper-time clocks and correlated sampled channels are covered
  by the research test suite.
- Additional convergence run:
  https://github.com/aalida0/QuirkSim/actions/runs/42835789643,
  three representative 2T `ng_endpoint` recoil events (seeds 10, 42, 197):
  max change in rmin under max_step=2 -> 0.5 ps is
  6.96e-9 um; under regulator 1e-7 -> 5e-8 cm is
  3.80e-7 um. This is a **spot check only**, not global convergence.

**Physical completeness: BLOCKED.** Even after incorporating the 2020
straight-string endpoint-force approximation:

1. The actual moving Nambu-Goto worldsheet PDE and its nonlocal
   propagation of transverse impulses are NOT evolved.
2. String energy and angular momentum are absent from the bookkeeping
   when the string is transversely moving; the endpoint force cannot
   replace their dynamics in all regimes.
3. Magnetic field is uniform, not a map; material ionization, scatters,
   infracolor/pion radiation and annihilation are absent.
4. Both endpoints are assumed charged at production, a conditional
   channel rather than an inclusive doublet/triplet production mixture.
5. The imposed CM boost is only one kinematic slice; real pp+jet event
   kinematics require a separate public EW-doublet/triplet UFO and MG5.
6. Weak decays are allowed to occur after repeated nominally
   core-overlapping encounters; without a dynamical contact survival
   model, predicted weak-decay fractions can be unphysical.
7. The significantly different **1 vs 3 um** zero-field recoil median
   in the two force prescriptions proves substantial model dependence
   even with numerically stable integrators.

## Physics verdict

**A magnetic field gives a physically motivated origin for pre-decay
nonradial motion**, unlike micron-scale initial endpoint offsets.
The recoil-induced geometric shift survives in both approximate force laws
in this kinematic slice. However, neither force prescription is a
validated full string simulation and the numerical distance predictions
differ by factors of several.

**PAPER: CONDITIONAL, not yet GO.** Next gate: independently verified
energy-momentum/ang-momentum-conserving dynamical string (worldsheet PDE or
well-characterized effective theory), including finite signal propagation,
then radiation and annihilation modeling with explicit uncertainties.
