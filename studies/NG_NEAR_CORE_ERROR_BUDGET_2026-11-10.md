# Near-contact endpoint-geometry error budget — 2026-11-10

**Status: ODE_AND_REGULATOR_STABILITY_PASS; SPATIAL_CORE_GEOMETRY_BLOCKED.**

This is a controlled numerical validation of the research
`EquidistributedNGString` discretization. It is **not** a physical
annihilation calculation or a complete 1000 ps worldsheet trajectory.

## Setup

- Initially pre-stretched Dirac triplet pair: neutral mass 225 GeV,
  Λ_IC=774 eV, relative beta=0.4, x CM boost=0.25.
- Uniform magnetic field B_z=2 Tesla.
- Zero **artificial** interior point masses; constrained equal-arclength
  semidiscrete Nambu–Goto action.
- Stop at the **same total spatial string arclength**:
  L_cut=2.0e-6 cm. Do not continue through this near-zero-length
  boundary because the physical continuation is unvalidated.
- Compare 5, 7, 11, 15 string vertices.
- For each vertex count:
  * nominal max_step=0.5 ps, rtol=2e-9, atol=1e-10,
    short-distance regularizer=1e-9 cm;
  * tight max_step=0.125 ps, rtol=2e-10, atol=1e-11;
  * half regularizer=5e-10 cm, otherwise baseline settings.
- All distances are in cm, not micrometres.
- Geometric reference length hbar*c/Λ_IC=2.54944e-8 cm is **not a
  quantum annihilation radius**. A diagnostic gate requires inter-mesh
  chord disagreement <0.1 times that reference, without suggesting the
  gate alone predicts any physical encounter probability.

## Independent numerical results

| Vertices | Endpoint chord, tight integrator (cm) | |tight − nominal| (cm) | |half-regulator − nominal| (cm) |
|---:|---:|---:|---:|
| 5 | 1.99806169e-6 | 2.78e-13 | 3.99e-15 |
| 7 | 1.99300800e-6 | 1.89e-12 | 3.59e-13 |
| 11 | 1.97071229e-6 | 1.46e-14 | 3.92e-13 |
| 15 | 1.99612611e-6 | 3.99e-15 | 3.09e-13 |

These are the values in the **combined** CI run; separately parallel
runs reproduce the results within numerical tolerance.

The **largest variation from the tighter integrator** is 1.893e-12 cm;
from the halved regulator, 3.924e-13 cm. Both are tiny compared with
the absolute **inter-mesh chord spread**:

```text
max(chord_n) - min(chord_n) = 2.73494047e-8 cm
hbar*c/Lambda_IC             = 2.54944419e-8 cm
inter-mesh spread/reference  = 1.07275950
```

Moreover, the vertex-count sequence **is not monotonic**:
n=11 has a notably smaller endpoint separation than n=5,7,15.
Tighter ODE tolerances do not resolve the mesh-level discrepancy.

## Diagnostic outcome

- **PASS:** Mesh-independent total arclength at the crossing; strict
  ODE tolerance and regulator stability at each tested mesh; energy
  conservation (max absolute energy drift <=6.82e-13 GeV).
- **FAIL/BLOCKED:** Core-reference chord convergence with mesh count;
  it misses the <=0.1 x reference scale requirement by more than a
  factor of 10.
- **BLOCKED:** Interpretation of endpoint closest approach as quantum
  annihilation; physical contact continuation; weak impulse on this
  zero-interior-mass gauge; all subsequent multi-return trajectories.

The reduction in global energy drift is **not** proof of geometric
convergence. The missing physics is not simply a smaller Runge–Kutta
step or smoother short-distance regulator. The proper next check is
interleaved mesh refinement and a gauge-/finite-element discretization
audit (including the midpoint quadrature and treatment of higher
transverse modes). Arbitrary remapping through zero length is not
authorized.

## Reproducibility

- [Combined stringent mesh/ODE/regulator run](https://github.com/aalida0/QuirkSim/actions/runs/42891059972)
- [Parallel per-mesh independent runs](https://github.com/aalida0/QuirkSim/actions/runs/42891915791)
- Implementation: `scripts/near_core_error_budget.py`.
- This research remains on `research/variational-worldsheet`; the
  stable development branch and 1000-event pilot data are untouched.

```bash
python -u scripts/near_core_error_budget.py \
  --nodes 5,7,11,15 --cutoff-fraction 0.00002
```

**Paper: still CONDITIONAL / NOT YET PHYSICS GO.**
