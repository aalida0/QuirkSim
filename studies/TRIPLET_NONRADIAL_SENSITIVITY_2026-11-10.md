# Triplet transverse-offset / initial-angular-momentum sensitivity

**Date:** 2026-11-10. **Status:** numerical sensitivity only; NOT a collider,
annihilation, radiation, or moving-string prediction.

## Setup

Standalone QuirkSim central straight-string **potential proxy**
(with smooth regulator), fully propagated charged-to-neutral pion, electron
and muon decays. Two initially charged endpoints with opposite longitudinal
momenta, beta=0.4, neutral quirk mass 225 GeV, Lambda_IC=774 eV, evolved
for 1000 ps. The physical recoil/no-recoil diagnostic share the same
proper-time clocks, decay-channel assignments and decay-rest-frame seeds,
even though lab-frame decay times can differ.

Transverse initial endpoint separation b=|x1-x2| (positive impact parameter)
implies **nonzero initial orbital angular momentum** L = b * pz.
This scan fixes beta, string tension, lifetime and random samples while
changing b. These b values are *not* realistic pp production-vertex offsets.

**Reference geometric core:** hbar*c/Lambda_IC ~ 0.000255 micrometres.
Encounters below this length scale are not annihilation predictions.

## Reproducible CI studies

- Exploratory 20 matched seeds per offset, b=0,0.01,0.1,1,3,10 micrometres:
  [successful run](https://github.com/aalida0/QuirkSim/actions/runs/42831898480).
- Follow-up 200 matched seeds per offset for b=0.1,1,3,10 micrometres:
  [successful run](https://github.com/aalida0/QuirkSim/actions/runs/42832418438).
- Raw JSON/CSV results are attached as workflow artifacts. The scan and
  regression tests are in `scripts/triplet_impact_scan.py` and
  `tests/test_triplet_impact_scan.py`.

**200-seed comparison**, in micrometres. Values are medians across events;
the paired increase is the median *event-level difference*, not the difference
of separate medians.

| initial b (um) | median no-recoil rmin (um) | median recoil rmin (um) | median paired increase (um) | fraction recoil rmin > no-recoil | fraction with geometric post-decay core contact (no/recoil) |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.100014 | 2.462289 | +2.362273 | 197/200 (98.5%) | 0/200 / 0/200 |
| 1 | 1.000140 | 2.717643 | +1.717436 | 172/200 (86.0%) | 0/200 / 0/200 |
| 3 | 3.000421 | 3.694916 | +0.694333 | 121/200 (60.5%) | 0/200 / 0/200 |
| 10 | 10.001404 | 10.001871 | +0.000345 | 100/200 (50.0%) | 0/200 / 0/200 |

**20-seed exploratory small-b comparison:**

| b (um) | no recoil rmin median (um) | recoil rmin median (um) | paired increase median (um) | geometric core contacts (no/recoil) |
|---:|---:|---:|---:|---:|
| 0 | ~0 | 2.831218 | +2.831218 | 20/20 / 0/20 |
| 0.01 | 0.010002 | 2.833590 | +2.823588 | 0/20 / 0/20 |

## What is and is not established

1. The pion/lepton weak-decay recoil creates a substantial trajectory
   perturbation relative to small pre-existing impact parameters, within
   the simplified central-string model. For b<=1 um this remains evident
   in the 200-seed medians and paired signs.
2. The geometric perturbation decreases as imposed initial angular momentum
   increases. At b=10 um the sign is balanced exactly 100/200, and median
   event-level increase is only 0.00035 um.
3. **The strong difference in geometric core-contact fractions occurs
   only for the ideal zero-impact-parameter reference.** At all tested
   nonzero offsets b>=0.01 um, both recoil and no-recoil trajectories have
   ZERO post-decay geometric core contacts in the sample.
4. This does not show a physical suppression of quirk annihilation. The
   no-recoil control is unphysical, annihilation is not simulated, and the
   simple core-radius threshold is not an annihilation rate.
5. A pointlike hard process produces both particles at the same location;
   in the pair CM, different opening angles alone do not create initial
   internal orbital angular momentum. Consequently the micron-range
   starting offsets in this test are **controlled stress inputs**, not
   a prediction of pre-existing orbital angular momentum at the LHC.
6. A realistic conclusion requires a physically motivated source of
   pre-decay angular momentum, proper moving Nambu--Goto string dynamics,
   radiation and prompt-contact survival, and collision kinematics.
   It is not sufficient to run more of the same b-offset samples.

## Implication for paper feasibility

**CONDITIONAL / NOT YET GO.** The simulation yields a reproducible recoil
geometry effect but does not yet establish a novel physical difference
in annihilation or radiation. These outputs should be considered
numerical/phenomenological prototypes, not figures for a submitted paper.
