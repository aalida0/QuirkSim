# Core-scale geometry versus mesh/regularizer — quirk string Step 1

**Date:** 2026-11-10  
**Status:** `PRECONTACT_MESH_REFINEMENT_PROGRESS / INFRARED_CONTACT_AND_CONTINUATION_NOT_VALIDATED`

This report concerns the **undecayed**, initially stretched Dirac-triplet
string (225 GeV neutral mass, Lambda_IC = 774 eV) in a uniform 2 T field
with a fixed beta_CM,x=0.25 and beta_rel=0.4. It is a solver benchmark,
NOT a collider-generated initial state or a weak-decay recoil result.

## Scientific question

Earlier N=5/7/11 meshes gave nearly identical *times* of first crossing
of a fixed total-string-arclength cutoff, but endpoint-chord differences
were comparable with the confinement scale. Since a curved/folded string
can have spatial arclength much greater than the chord between its
endpoints, a converged contact time is **not evidence of a converged
core geometry**.

We separately vary three numerical controls:

1. Spatial node count N, with zero artificial interior masses and
   constrained equal-spatial-arclength gauge;
2. The spatial string-segment length regularizer a;
3. DOP853 maximum integration step.

The reference length is hbar*c/Lambda_IC =
`2.549444191783476e-8 cm`; **it is a geometric diagnostic, NOT a
quantum annihilation radius**. The initial string length is 0.1 cm.

## Controlled fixed-arclength comparison

At a common TOTAL spatial string arclength `Lcut=2e-6 cm` (about
78.45 reference confinement lengths), with
`a=1e-10 cm, max_step=0.5 ps`:

| N | Endpoint chord (cm) | Lcut minus chord (cm) | Maximum adjacent segment bend (rad) |
|---:|---:|---:|---:|
| 11 | 1.9707118243e-6 | 2.92881757e-8 | 0.09551 |
| 15 | 1.9961256922e-6 | 3.87430831e-9 | 0.05527 |
| 21 | 1.9993889794e-6 | 6.11020913e-10 | 0.02186 |
| 25 | 1.9996246565e-6 | 3.75341478e-10 | 0.01724 |

- The 11→15 difference is `2.5414e-8 cm` = **0.997** reference lengths.
- The 15→21 difference is `3.2633e-9 cm` = **0.128** reference lengths.
- The 21→25 difference is `2.3568e-10 cm` = **0.00924** reference lengths.

The last three mesh differences shrink; **this is a provisional fine-grid
chord-convergence PASS at the stated cutoff only**. It is not a
physical annihilation statement.

### Independent integrator and regularizer scans

At N=21 with `a=1e-10 cm`, changing max_step from 0.5 to 0.2 ps
changes the chord by `1.33e-15 cm`, far below the 1/Lambda scale.

At fixed N=11 or N=15, reducing the regularizer tenfold,
`1e-9 -> 1e-10 cm`, changes the chord by less than
`5e-13 cm`. The original large chord discrepancy is therefore
predominantly **spatial under-resolution**, not time-step or regularizer
error in this comparison.

The original 1e-9 cm regularizer plus its conservative per-element
mesh guard blocked N=21 and N=31 prematurely; one must NOT compare
their stopped states at *different total string lengths* as if they
belonged to the same physical encounter. With 1e-10 cm regularization,
N=21 and N=25 reach the common physical total-length event cleanly.

## Much closer to the infracolor scale: independent gate

The next target is total spatial arclength `Lcut=1e-7 cm`,
about 3.92 times hbar*c/Lambda_IC. At `a=1e-12 cm`,
the original velocity-coordinate constrained NG integrator is
**BLOCKED: superluminal intermediate string vertex** (N=11 and N=15).
This is a numerical coordinate failure, **not physical superluminality**.

A separate solver combines the zero-interior-mass action with
bounded spatial rapidity variables w=gamma*v (so |v|<c) and
retains the same arclength, mesh and conservation events.
This formulation is also an experimental pre-contact method;
its success or failure at an Lcut says nothing about the physical
continuation through a string collision. Any near-core outcomes
must be independently checked against mesh and integration
refinements before inclusion in a phenomenology paper.

## Near-infracolor gate at ~4 reference confinement lengths

A much closer test uses total-arclength cutoff
`Lcut=1e-7 cm` (3.92 times hbar*c/Lambda_IC), spatial regulator
`a=1e-12 cm`, 2T and the same boosted pre-stretched initial slice.

**Velocity-coordinate constrained solver:** the N=11 and N=15
calculations are **BLOCKED by a superluminal intermediate numerical
vertex**. This is an invalid integration stage, not evidence of a
physical faster-than-light trajectory.

**Bounded-speed rapidity constrained solver:**

| Method | Last valid arclength (cm) | Last valid time (ps) | Status |
|---|---:|---:|---|
| DOP853 | 2.823751e-7 | 458.6120574493 | Step size below float resolution |
| Radau | 2.826045e-7 | 458.6120574458 | Step size below float resolution |
| BDF | 5.487967e-7 | 458.6120419879 | Step size below float resolution |
| LSODA | not certified | not certified | velocity domain error |

No solver has demonstrated the prescribed `Lcut=1e-7 cm`
trajectory. All these runs **stop before** the requested cutoff.

**Disambiguation of worldsheet versus gauge speeds.** At the last valid
DOP853 state:

- Maximum **individual numerical vertex speed**: `0.999993799c`.
- Maximum **segment transverse speed**:
  `0.932873923c`.
- Minimum target-time-gauge segment lapse squared
  `1-|v_perp|² = 0.129746243`.
- Maximum absolute segment-parallel speed `0.442346313c`.

Since the transverse segment lapse remains **positive and not near
zero**, there is NO demonstrated physical Nambu-Goto cusp at this
failure. The problem is consistent with a singular or extremely
stiff spatial parametrization. This makes a coordinate/gauge
reformulation better justified than applying an arbitrary
velocity or curvature cap.

The endpoint chord is `1.62643869e-7 cm` while the total spatial
arclength is `2.82375103e-7 cm`; the string has substantially
nonstraight geometry on this last valid slice. Neither number
is a physical annihilation probability.

**Evidence:** [bounded-speed constrained NG](https://github.com/aalida0/QuirkSim/actions/runs/42859333149)
and [Radau/BDF/LSODA cross-check](https://github.com/aalida0/QuirkSim/actions/runs/42858633861).

## Physical origin of possible folds

Sonnenschein, Weissman and Yankielowicz,
[The scattering amplitude of stringy hadrons: strings with opposite
charges on their endpoints](https://link.springer.com/article/10.1007/JHEP07%282020%29156)
(JHEP07(2020)156), show that oppositely electrically charged string
endpoints in a constant magnetic field can admit rotating folded
open-string solutions. Their particular classical solutions and mass
assumptions are **not** directly this massive-quirk event, but they
show that folds are potentially physical. An equal-arclength gauge
cannot be blindly continued through a gauge-degenerate pinch.

## Physics verdict

**PASS for this one pre-stretched, no-decay 2 T slice:** zero artificial
interior point masses, numerical energy conservation, nearly identical
event-hitting times, short-time rapidity/velocity equivalence, and a
decreasing fine-grid chord discrepancy at Lcut=2e-6 cm.

**BLOCKED for paper use:** validated continuum worldsheet and a proper
massive-endpoint crossing/continuation, core-scale chord convergence
at physically relevant shorter distances, weak decay on the
new zero-mass constrained gauge, radiation/annihilation and survival,
and correct hard-process initial conditions.

A successful CI job only means the diagnostic executed. It does not
override these scientific blockers.

## Reproduce

```bash
git switch research/variational-worldsheet
python -m pip install -e '.[test]'

# Fixed 2e-6 cm spatial arclength at N=21 and N=25
python scripts/run_near_contact_shape.py --nodes 21 \
  --regulator-cm 1e-10 --pinch-fraction 0.00002 \
  --output results/near_contact_n21.json
python scripts/run_near_contact_shape.py --nodes 25 \
  --regulator-cm 1e-10 --pinch-fraction 0.00002 \
  --output results/near_contact_n25.json

# Smaller core-scale cutoff, bounded-speed coordinate formulation
python scripts/rapidity_constrained_contact.py --nodes 11 \
  --regulator-cm 1e-12 --pinch-fraction 0.000001 \
  --output results/rapidity_n11.json
```

**CI evidence:**
- [Short-distance regulator sensitivity (N11/N15/N21)](https://github.com/aalida0/QuirkSim/actions/runs/42856016395)
- [N25 fine-grid chord run](https://github.com/aalida0/QuirkSim/actions/runs/42856564224)
- [Original/coarse-grid test and explicitly blocked mesh guard](https://github.com/aalida0/QuirkSim/actions/runs/42855532002)
