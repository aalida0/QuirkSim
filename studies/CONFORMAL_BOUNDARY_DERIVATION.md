# Temporal conformal massive-boundary prototype

Date: 2026-11-10. This implements the first portion of the conformal-solver
specification. It does **not** complete the magnetized 1000 ps worldsheet,
contact prescription, weak recoil, or survival calculation.

## Source and independent specialization

The moving-boundary action and flux in Avgoustidis, Pourtsidou and
Sakellariadou, *Phys. Rev. D* **91**, 025022 (2015),
[DOI](https://doi.org/10.1103/PhysRevD.91.025022), section III,
equations (20)–(35), supply the starting point. The paper treats a massive
junction. The one-leg endpoint specialization, orientation signs, units,
and incoming-ray elimination below are derived here. The implementation
is not a reproduction of a published numerical endpoint code.

## Coordinates and endpoint flux

Use laboratory time t in ps, an energy coordinate s in ps, and y=x/c
in ps. Write y=[a(t+s)+b(t-s)]/2. Define A=a', B=b', both unit vectors.
Then y_t=(A+B)/2 and y_s=(A-B)/2; the wave equation and both temporal
conformal constraints hold on every smooth bulk patch. Our b argument
is the negative of the paper's convention; its derivative changes sign.

The string has endpoint labels s_L(t)<s_R(t). For
f=Lambda_IC^2/hbar in GeV/ps, variation of the action gives

```
dp_L/dt = +f (y_s + ds_L/dt y_t)
dp_R/dt = -f (y_s + ds_R/dt y_t)
dE_L/dt = +f ds_L/dt
dE_R/dt = -f ds_R/dt
```

Here v_i=p_i/sqrt(m_i^2+|p_i|^2)=y_t+ds_i/dt y_s is the physical
endpoint beta. These are ordinary endpoint momenta, not canonical
momenta attached to finite-element labels. In this first implementation
external electromagnetic forces are absent.

At the left end take incoming K=A, outgoing O=B, and q=ds_L/dt.
At the right end take incoming K=B, outgoing O=A, and q=-ds_R/dt.
For both ends boundary kinematics and the outgoing null constraint give

```
2v = (1+q)K + (1-q)O
q  = (v.K - |v|^2)/(1-v.K)
O  = [2v-(1+q)K]/(1-q)
dp/dt = f [(1+q)K-v]
dE/dt = f q
```

The last identity also follows by contracting the momentum equation
with v. It fixes the energy-transfer sign without an assumed damping
term. In the radial sector K is the inward unit tangent, the force is
exactly fK, and the endpoint labels follow longitudinal endpoint motion.

## Independent charge and geometry diagnostics

At a common laboratory time, the history integrals give

```
E_string = f (s_R-s_L)
P_string = f/2 [integral_(t+s_L)^(t+s_R) A du
                 + integral_(t-s_R)^(t-s_L) B dv]
x_R-x_L = c/2 [integral A du - integral B dv]
```

In zero external field the sum of endpoint and string energy is conserved:
the string energy derivative is -f(q_L+q_R), opposite the endpoint sum.
Momentum conservation and chord closure are evaluated using history
integrals, independently of the endpoint position ODE quadrature.
Angular momentum and boost-charge diagnostics are not implemented yet.

## Numerical implementation and limitations

`conformal_boundary.py` contains the algebraic massive-boundary map.
`conformal_characteristics.py` propagates two characteristic histories
and advances boundary positions, momenta, labels and proper times with
RK4. Each step is below 0.2 times the current label interval, ensuring
all incoming rays are already available from initial data or earlier
endpoint emission. Missing history raises an error instead of extrapolating.

Directions between stored samples follow great-circle arcs. Their
primitives are integrated analytically. Analytically unit outgoing rays
are normalized only after rejecting deviations greater than 1e-9, to
avoid accumulated roundoff under repeated reflection. This construction
enforces the bulk null constraints; tiny constraint residuals alone do
not validate the endpoint trajectory. Initial-data compatibility, independent
analytic trajectories, conserved charges, and independent time/initial-grid
refinement are additional checks.

The generic solver stops when its energy-coordinate interval approaches
zero. It neither sorts endpoint positions nor installs a reflection or
annihilation rule. A signed radial adapter in `conformal_radial.py`
separately reverses the ideal string orientation at each zero-length
crossing, keeps the endpoint momenta continuous, and records the chosen
prescription. It uses the boundary flux and numerical event integration;
it never calls `ExactRadialNG`. The latter is an independent analytic
comparison in the validation script. This adapter is a symmetry reduction
with extra continuation data, not evidence for general 3D continuation.

## Independent nonlinear rotating benchmark

Besides radial motion, use

```
y(t,s) = sin(omega s)/omega (cos(omega t), sin(omega t), 0)
s_end = +/- asin(beta)/omega
omega = f (1-beta^2)/(m beta)
```

This solves the bulk wave equation, both constraints, and the massive
boundary equation: m gamma beta omega=f sqrt(1-beta^2).
Initial A=(cos(omega s),sin(omega s),0) and
B=(-cos(omega s),sin(omega s),0) generate the exact solution.
The endpoints have radius c beta/omega, opposite tangential momenta,
and proper time t sqrt(1-beta^2). The benchmark has transverse string
motion and massive dynamical endpoints; it is distinct from the radial
constant-force check. It remains a special symmetric solution.

## Reproduction and release decision

```
python -m pip install -e '.[test]'
python -m pytest -q
python scripts/validate_conformal_boundary.py
```

The validation script predeclares numerical tolerances and emits JSON
with separate radial, rotation, pulse-time and initial-grid checks.
See `CONFORMAL_BOUNDARY_RESULT_2026-11-10.md` for witnessed values.

**Release:** research prototype only. Zero-field radial and rotating
analytic checks can be reported with their exact scopes. General
transverse convergence, all Poincare charges, electromagnetic boundaries,
the selected finite-element cross-comparison, curved near-contact
continuation, physical weak recoil and radiation/annihilation remain
unvalidated. No new geometric-contact fraction or collider sensitivity
may be inferred from these checks. Keep this on the research branch.
