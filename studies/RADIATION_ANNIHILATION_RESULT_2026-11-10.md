# Causal evolution through radiation/annihilation: implementation and physical limits

Date: 2026-11-10. Research branch only; no merge to the stable simulator or main.

**The conditional radiation/annihilation machinery is implemented. Absolute
physical survival and recoil-induced annihilation suppression are still UNKNOWN.**
The remaining inputs are quantum/nonperturbative matching data, not missing
ODE tolerances. Successful code execution cannot supply those data.

## What changed

The conformal characteristic solver now includes charged massive endpoints in
a uniform magnetic field, field impulse/torque/boost accounting, and classical
orbital Poincare diagnostics. Bulk string charge is zero. Weak pion, electron
and muon transitions conserve four-momentum locally and update the endpoint
mass and electric charge. The outgoing characteristic has a one-sided jump;
the corresponding signal and its reflections travel causally. The incoming
string and the bulk do not receive a nonlocal instantaneous momentum update.

Close returns exposed an error that energy conservation alone missed.
Independent causal-fraction and ray-angle refinement now controls it.
The solver also locates an explicit short-distance matching surface rather
than equating a radius crossing with annihilation. The surface censors the
classical calculation unless an explicit quantum rate/population card is supplied.

The new radiation/annihilation module provides:

- Smooth single-charge Lienard power in an explicitly derived EM convention.
  The 3D photon ledger is a diagnostic without radiation reaction; its entries
  are not counted again as actual emitted energy in the mechanical ledger.
- Same-flavor, color-neutral, unpolarized Dirac IC-singlet hard coefficients
  for IC gluons and QED channels, with threshold factors and explicit coverage.
  They are **partial coefficients**, not the complete EW doublet/triplet total.
- A radial zero-field orbit average with supplied glueball energy/efficiency,
  energy-loss feedback, supplied s-wave overlap, and competing first weak decay
  and annihilation probabilities. Surviving until a later weak decay is accounted
  for explicitly; killed paths do not continue contributing probability.
- A finite-level quantum master equation with radiative transitions, annihilation
  and first-weak sinks. Probability and mean energy have separate conservative
  ledgers in that fixed-center model. This is not an event generator with
  exact discrete radiation recoil. Level energies, populations, allowed mass
  thresholds, selection rules and transition rates require physical matching.
- An explicit classical-to-quantum handoff that matches invariant energy and
  uses the correct proper time. It assumes a fixed boost and isotropic CM mean
  radiation after matching. No benchmark receives an invented rate card.

## Witnessed validation and evidence

The complete local suite passes **135 tests**. Analytic checks include an
angular Larmor integral, the constant longitudinal-force radiation limit,
parapositronium normalization, explicit SU(N) color traces, a conservative
two-stage radiation-to-annihilation cascade, exact competing exponentials,
independent cooling-time quadrature, charge conservation and causal weak-wave
reflection. These certify the stated calculations and implementations, not an
unprovided physical density matrix.

| Claim | Evidence | Decision and scope |
|---|---|---|
| 2 T magnetic endpoint evolution | FEM Gauss-2 with 7/11/15 vertices versus characteristics at identical lab times 20/100/200 ps | PASS for these pre-contact times; largest 15-node endpoint difference 2.495e-9 cm |
| Closer approach in the new coordinate system | Causal-fraction/ray-angle refinement to an 8 hbar*c/Lambda endpoint matching radius | PASS for selected matching geometry; finest chord-closure error 5.04e-10 cm, last hitting-time shift 3.71e-8 ps |
| Orbital angular/boost diagnostics | Exact rotating reference and field-source ledgers | Implemented and checked classically; quantum partial waves remain undetermined |
| Local weak decay and reflected null waves | All pion/e/muon channels on both endpoints; boosted early-decay regression | PASS for local four-momentum/charge matching and selected causal propagation |
| Generic full classical collision/reconnection | No quantum return map or reconnection rule supplied | UNKNOWN; matching stops instead of assigning survival |
| Smooth EM radiation | Explicit rationalized power convention and analytic tests | PASS for diagnostic power; pair interference and impulsive vertex radiation absent |
| Glueball radiation | Supplied epsilon/energy prescription and conservative orbit-average feedback | CONDITIONAL; calibrated differential emission/recoil law UNKNOWN |
| Hard QED/IC annihilation coefficients | Table-II source match, spin projection and independent color algebra | PASS for stated partial kernels; complete EW/flavor matrix UNKNOWN |
| Quantum radiation/annihilation evolution | Analytic finite-level cascade, probability/energy ledgers and proper-time handoff | PASS as a conditional solver; physical rate/population inputs UNKNOWN |

Primary numerical records:

- [Magnetic matching refinement](data/conformal_magnetic_matching_2026-11-10.json)
- [Independent FEM comparison](data/conformal_fem_comparison_2026-11-10.json)
- [Conditional radiation/annihilation study](data/radiation_annihilation_assumption_study_2026-11-10.json)

The requested 1000 ps stress pilots have independent geometry/energy/momentum
acceptance checks in `validate_conformal_evolution.py`. A time limit reached
with a failed check remains a failed numerical result. An initially completed
seed-25 recoil run had max chord-closure error 5.31e-9 cm, above the fixed
0.1 hbar*c/Lambda tolerance; its tiny energy residual did not rescue it.
Refinement attempts are preserved in the JSON output. See the saved pilot
summary and CI artifacts for the final bounded runs, including failed attempts.

The seed-25 strict 1000 ps recoil pilot also passes independently on GitHub
Actions (code `e467aa4`, run `42911265905`). The seed-1 CI job was canceled
at its original 20-minute execution budget, after reaching 907.30 ps with
chord-closure error 5.13e-10 cm. This is not recorded as a numerical pass.
Full 1000 ps stress jobs now run via `workflow_dispatch` with a 45-minute
budget; routine pushes run the component and source/formula checks. Local
long-run results and their explicit statuses are preserved separately.

The event locator was also optimized to reuse identical RK trial states and
the common first stage. A boosted 30 ps seed-25 before/after check has
identical serialized trajectories, decays, reflections and ledgers, while
measured runtime decreased from 4.19 s to 1.19 s on the same host. This is
an engineering equivalence check, not a change of the physical model or
acceptance tolerance; the [record is saved](data/conformal_cache_equivalence_2026-11-10.json).

## Conditional timescale result

The assumption study uses pre-stretched initial strings of rest length 0.1 cm
and endpoint beta=0.4. These are not pointlike pp production initial states and
are not the earlier straight central-string ensemble. NIC=2,3 and alpha_IC=0.04
are supplied assumptions. The glueball efficiency is 0, 0.1 or 1 per full signed
period; its phenomenological emission energy is Lambda_tension. The assumed
s-wave overlap is 0, 0.1 or 1. No overlap is fitted to a classical miss distance.

The reference QED coefficients use the module's supplied alpha_EM=alpha(0),
not an inferred hard-scale running coupling. The orbit average combines an
exact radial period with a leading excited-state WKB annihilation relation;
relativistic corrections to that relation and spin-resolved evolution are not
certified. The quoted digits describe the conditional calculation, not its
physical precision.

| Benchmark | Weak lifetime of one endpoint | Shortest linear cooling time in tested prescriptions | Largest partial-model probability of annihilation before first weak decay |
|---|---:|---:|---:|
| Triplet, neutral mass 225 GeV, Lambda_tension 774 eV | 258.52 ps | 1.685e10 ps | 3.03e-4 |
| Doublet, neutral mass 150 GeV, Lambda_tension 2154 eV | 45.69 ps | 8.012e8 ps | 5.54e-4 |

**Interpretation:** within these declared prescriptions the charged sector
usually undergoes weak decay long before radiative cooling to a small bound
state. The partial-model probabilities are not upper bounds on the physical
EW-multiplet annihilation rate: omitted channels, interference, spin/flavor
matching and the initial density matrix can change them.

For an orbit conditioned to remain in that charged sector for 1000 ps, the
largest tested losses are only 1.73e-6 GeV (triplet) and 4.17e-5 GeV (doublet).
However, assigning such emitted momentum to a heavy endpoint gives free-flight
recoil displacement scales of about 9 and 908 confinement reference lengths.
These are materiality scales, not bounds on nonlinear string motion. **Small
energy loss does not justify ignoring radiation when predicting a tiny
annihilation-sensitive miss distance.** The orbit-averaged model is therefore
not promoted into a validated 3D radiation-backreaction prescription.

## Why a unique physical annihilation result is not finished

1. A differential glueball emission law is missing: probability, mass spectrum,
   energy, direction, spin/angular momentum and the causal string response.
2. The quantum spin/flavor/partial-wave density matrix at the matching surface
   is missing. Classical total orbital angular momentum is not an s-wave
   occupation. The numerical orbital-charge uncertainty also greatly exceeds
   one hbar near the matching surface even when macroscopic geometry converges.
3. The complete EW short-distance annihilation matrix is missing. After one
   charged-to-neutral transition the flavor sector changes; the same-flavor
   IC-gluon coefficient cannot be applied unchanged to that state.
4. Physical quantum transition/capture/return rates are missing. The implemented
   finite-level equation requires those rates; it does not derive them from NG.

The inspected published prescriptions do not close these inputs for this model. The
[source audit](RADIATION_ANNIHILATION_SOURCE_AUDIT_2026-11-10.md) records the exact
formulas, normalization checks, period convention and model mismatches.
Continuing a classical trajectory through a radius and sampling an arbitrary
Bernoulli annihilation flag would add an assumption, not finish the physics.

## Reproduction and release

```
python -m pip install -e '.[test]'
python -m pytest -q
python scripts/validate_radiation_annihilation.py
python scripts/validate_conformal_magnetic.py
python scripts/validate_conformal_evolution.py --mode stable
python scripts/validate_conformal_evolution.py --mode recoil --seed 25
```

The last commands are bounded numerical stress checks and can fail their
fixed geometry gate. They preserve the final state, validation and refinement
attempts. Unmatched quantum survival is always `null`/`UNKNOWN`, never zero or
one by default. Source code, selected classical benchmarks and **explicitly
conditional** assumption studies are reviewable. Absolute survival,
recoil-induced annihilation suppression, physical branching fractions and
collider reach are **not released as predictions**.
