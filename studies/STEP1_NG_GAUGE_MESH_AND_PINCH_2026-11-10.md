# Step 1: gauge-constrained Nambu--Goto string / pre-pinch validation

**Date:** 2026-11-10
**Branch:** `research/variational-worldsheet`
**Status:** `PRE_PINCH_CONVERGENCE_VALIDATED_IN_SELECTED_BENCHMARKS`;
`PINCH_CONTINUATION_BLOCKED`; `LARGE_BOOST_STIFFNESS_BLOCKED`;
`FULL_PHYSICAL_WORLD_SHEET_NOT_VALIDATED`.

## Scope and motivation

This independent open-source project is testing electroweak triplet weak
recoil in macroscopic infracolor strings. The production project and
detector simulation are *outside* this numerical stage.

The earlier fixed-material-vertex spatial finite-element model had
artificial positive masses on internal string vertices. Near its first
return those vertices became extremely unevenly distributed; rapidity
variables removed unphysical superluminal velocity stages but did not
remove a collapsing mesh. Energy conservation alone did NOT validate
physical trajectories through the collapse.

Current triplet numerical slice: m_neutral=225 GeV, Lambda_IC=774 eV,
charged endpoints, beta_rel=0.4; uniform B_z=0 or 2T; selected
beta_CM,x=0.25 unless varied; string initially PRE-STRETCHED 0.1cm.
The boosted starting slice is not an exact Lorentz transformation of
a worldsheet history. No weak decays are included in the new
equal-arclength solver tests yet.

## 1. Independent exact massive-endpoint radial reference

Module: `quirksim/radial_ng_reference.py`.
CI: https://github.com/aalida0/QuirkSim/actions/runs/42844612951

The one-dimensional strictly collinear radial sector has a
piecewise-exact solution of the reduced NG action, with string energy
`T*|x_R-x_L|`, relativistic endpoint kinetic terms, and signed
crossing continuation. It has NO transverse string dynamics or
infracolor quantum annihilation.

With the benchmark above and initial separation 0.1 cm:

| Observable | Exact radial result |
|---|---:|
| First crossing | 220.03789114 ps |
| Subsequent crossing interval | 224.14230831 ps |
| Number of signed crossings within 1000ps | 4 |
| Maximum endpoint separation | 1.45110286 cm |
| Sampled relative energy drift | 0 at reported precision |

Analytic derivatives of position, momentum, and proper time are tested.
Relabeling interior radial positions by different monotone fractions
while retaining the same physical endpoints and no artificial masses
leaves the exact string energy and total momentum invariant.

This validates a **restricted 1+1D sector**, not a curved worldsheet.

## 2. Location of original fixed-mesh breakdown

In 0T without boost, the original five-vertex mesh first violates
its element-length-ratio guard at **217.96061652 ps**, still
**2.0773 ps before** the exact first crossing. The old seven-vertex
mesh violates it at **217.28840003 ps**. Endpoint separations are
0.0510 and 0.0673 cm, respectively: far from an actual contact.
The endpoint separation matches the independent exact reference at
those times to ~1e-14 cm, but the INTERIOR coordinate mesh is failing.

At B=2T with beta_CM,x=0.25 the old 5/7-vertex mesh imbalance
occurs around 224.3-224.9 ps. This is not proof of a physical
annihilation, radiation, or worldsheet pinch.

Run: https://github.com/aalida0/QuirkSim/actions/runs/42844612951

## 3. Gauge-constrained dynamical string WITHOUT internal gauge masses

Module: `quirksim/equidistributed_ng.py`.

The new semidiscrete NG solver sets the unphysical masses at all
interior material vertices **to zero**. To remove the resulting
tangential null directions in the canonical mass matrix, it imposes
the explicit gauge constraints

```text
g_j = |x[j+1]-x[j]| - |x[j]-x[j-1]| = 0,
d(g_j)/dt = 0,
d²(g_j)/dt² = J_j*(dv/dt) + b_j = 0.
```

It solves the variational saddle system

```text
[ M J^T ] [dv/dt] = [ physical force - convective correction ]
[ J  0  ] [lambda]  = [ -b                                      ].
```

Here `M=dP/dv` is built from the finite-element Nambu--Goto action,
including its internal string momentum. The gauge forces perform no
work when the constraints hold.

This is a **spatial midpoint finite-element NG approximation with
gauge constraints**, not an independently validated conformal-gauge
continuum worldsheet, and it cannot automatically traverse an
all-endpoints zero-length pinch.

### Initial 260ps B=2T validation

Both n=5 and n=7 reach **260ps** while maintaining an equal-length
mesh and energy to numerical precision:

| Quantity | n=5 | n=7 |
|---|---:|---:|
| Endpoint separation at 260ps | 0.67055187502 cm | 0.67055187505 cm |
| Largest element ratio at 260ps | ~1.0000000001 | ~1.0000000002 |
| Absolute energy drift | ~4e-13 GeV | ~0 GeV |
| Artificial internal gauge masses | **zero** | **zero** |

At 200ps a separate n=5/7/11 endpoint-position comparison shows
maximum difference versus n=11 of 9.09e-8 cm; maximal energy drift
is 2.84e-13 GeV.

Runs:
- https://github.com/aalida0/QuirkSim/actions/runs/42847499209
- https://github.com/aalida0/QuirkSim/actions/runs/42848882971

## 4. Crucial correction: a COMMON TOTAL arclength stop

A stopping threshold imposed per mesh element is NOT invariant under
changing the number of vertices. It artificially stopped 11 vertices
before n=5/7. The correct controlled numerical comparison uses one
common **total string polyline arclength** cutoff, independent of n,
plus separate local element-quality safety stops.

The comparator
`scripts/check_total_arclength_convergence.py` and unit tests
`tests/test_total_arclength_gate.py` explicitly enforce this.

At fixed B=2T, beta_CM,x=0.25, the n=5,7,11 stopping times are:

| Common total string arclength cutoff | Approx. stop time (ps) |
|---|---:|
| 4e-5 cm | 458.61036 |
| 4e-6 cm | 458.61189 |
| 1e-6 cm | 458.61202 |

The spread in stop time across n=5,7,11 is only O(1e-6ps)
at the smaller cutoffs, and the combined string energy drift
remains below O(1e-12GeV). The three resolutions are approaching
a VERY SHORT configuration near t≈458.612ps.

**Critical qualification:** The stop itself is an engineering event;
it does NOT establish exact zero string length, self-intersection,
annihilation, or a unique continuation. Even at 1e-6cm, the
stopping arclength is about 39 times the reference length
`hbar*c/Lambda_IC` and cannot be used as an annihilation threshold.
The post-cutoff trajectory is **not** propagated.

Independent strict two-cutoff (1e-4 and 4e-5 cm) run:
https://github.com/aalida0/QuirkSim/actions/runs/42850022362
(all numerical convergence gates PASS). Smaller cutoff runs:
https://github.com/aalida0/QuirkSim/actions/runs/42849334543

## 5. Higher CM boost / subluminal-coordinate robustness

`scripts/equidistributed_rapidity_stress.py` combines the zero-mass
constrained gauge with the rapidity map
`v=w/sqrt(1+w.w)`. At B=2T, total arclength cutoff=1e-6cm,
and requested horizon=1000ps:

| beta_CM,x | n=5 | n=7 |
|---|---|---|
| 0.25 | stop at ~458.612ps, good energy | same stop within ~2e-6ps |
| 0.45 | stop at ~496.911ps, good energy | same stop within ~5e-6ps |
| 0.65 | stops at ~583.333ps, energy drift ~1e-12GeV | **BLOCKED**, integrator step underflow |
| 0.80 | **BLOCKED**, integrator step underflow | **BLOCKED**, integrator step underflow |

Neither of the high-boost numerical failures is evidence for a
superluminal physical string: the rapidity map keeps vertex
speeds below c. The failures mean the finite-element/KKT
evolution is too stiff or singular to be trusted there.

Run:
https://github.com/aalida0/QuirkSim/actions/runs/42850434436

## 6. Scientific verdict and required next steps

**Achievement:** The previous spurious interior-mesh collapse is
separated from physical endpoint motion, the radial 1000ps sector has
an independent exact reference, and a new zero-artificial-mass
variational gauge solver provides strongly convergent **pre-pinch**
3D numerical trajectories for selected 2T kinematic slices.

**NOT ACHIEVED:**
- A continuum-validated conformal gauge/worldsheet evolution through
  a genuine zero-length pinch, collision or kink.
- Guaranteed causal transport and stable reparameterization of folds.
- A convergent worldsheet solver at all CM boosts.
- A working weak-decay canonical impulse coupled to the NEW zero-mass
  constrained gauge: the previous weak impulse was tested for a
  different finite-element model with artificial internal masses and
  cannot simply be re-used.
- Exact Lorentz-transformed physical pp initial conditions.
- Self-consistent infracolor radiation or quantum annihilation rates.

**Step 1 status: significant progress, NOT COMPLETE.**
The next correct mathematical milestone is a massive-endpoint
conformal/characteristic formulation or a carefully derived
collision/continuation boundary prescription, independently
benchmarked against the exact radial reference. A physically
unjustified remesh at contact would only hide the problem.

**Paper status: CONDITIONAL.** Neither a 1000ps 3D NG evolution nor
an annihilation-suppression claim has passed physics validation.
