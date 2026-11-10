# Next Step-1 gate: independent conformal-gauge massive-endpoint worldsheet

**Scientific state (2026-11-10):** the pre-contact Gauss-quadrature
constrained finite-element worldsheet PASSES selected triplet near-core
geometry refinement. **No physical post-contact Nambu-Goto worldsheet
continuation has been established**, and the 1000 ps 3D study is
NOT complete.

## Why the current spatial gauge cannot be blindly continued

In target-time gauge with equal spatial arclength, a very small total
length makes the tangential constraint matrix ill-conditioned.
The current code intentionally stops at that boundary. Neither a
successful ODE solve nor conserved energy validates extrapolation
through a zero-length or self-intersecting worldsheet.
Confinement scale 1/Lambda_IC is a reference length; it does not
supply the missing quantum annihilation boundary condition.

A mathematically different parameterization is required as an
**independent cross-check**, not an ad-hoc relabeling of collapsed
vertices. The natural candidate is a 1+1-dimensional conformal
worldsheet coordinate system with the wave equation

```text
(∂_τ² - ∂_σ²) X^μ = 0
(∂_τ X ± ∂_σ X)^2 = 0
```

and *massive dynamical endpoint boundary conditions* derived from
the Nambu-Goto + endpoint-particle action. The endpoint
equations balance incoming string momentum flux and, in the
magnetized quirk problem, an external Lorentz four-force.
The string can carry its own energy and angular momentum.

The historical boundary-value literature specifically studies
relativistic strings with finite endpoint masses and analytic
extension of initial data by endpoint ODEs. This is a
**research lead to be independently reproduced**, not proof that
arbitrary macroscopic quirk self-intersections have a unique physical
reconnection/annihilation prescription.

## Literature to derive and audit before coding

1. B. M. Barbashov and G. S. Sharov, *Initial-boundary problem
   for the relativistic string with massive ends*, Theor. Math.
   Phys. **101** (1994) 1332–1345,
   DOI 10.1007/BF01018281.
   https://www.mathnet.ru/php/archive.phtml?jrnid=tmf&option_lang=eng&paperid=1683&wshow=paper
   Abstract: general initial-position/velocity boundary problem,
   reparameterization invariance; finite-mass extension reduces
   to an ODE boundary system.
2. G. S. Sharov, *Solution of the initial boundary value problem for
   the relativistic string with massive ends*, Comput. Math. Math.
   Phys. **37** (1997) 590–601.
3. B. M. Barbashov and V. V. Nesterenko,
   *Relativistic string with massive ends*, Theor. Math. Phys.
   **31** (1977) 465–470. An orthogonal-gauge massive endpoint
   construction and Fourier-mode conditions.
4. Open-access manuscript on quirk strings:
   *The quirk trajectory*, EPJC **80** (2020) 651,
   https://link.springer.com/article/10.1140/epjc/s10052-020-8209-y
   This provides an approximate straight-string endpoint-force
   benchmark but is NOT a full moving-worldsheet solution.

## Proposed independent validation contract

Stage A: **zero field, radial symmetry**. Reproduce
`ExactRadialNG` with identical endpoint initial data through
four signed crossings in 1000 ps, including endpoint proper time,
energy and momentum, without imposing an annihilation hazard.
The signed crossing is a *declared classical continuation convention*
rather than a quantum survival probability.

Stage B: **small transverse mode, no magnetic field**.
Reproduce the linearized worldsheet propagation at speed c,
conformally constrained mode equations, massive endpoint response
and conserved string+endpoint Poincare charges. Refine spatial
frequency and time resolution independently.

Stage C: **2 T boosted triplet, no weak decays**.
Compare conformal-gauge endpoint trajectories and string
energy-momentum to the selected `GaussQuadratureNGString`
13/17/21 vertex pre-contact data. Require changes in endpoint
separation below a predeclared tolerance relative to
hbar*c/Lambda_IC; verify gauge constraints at every recorded step.

Stage D: **contact boundary**. Classify whether a continuation
is regular in the conformal formulation, genuinely singular
or dependent on extra self-intersection/reconnection data.
Do not silently relabel or reflect endpoints. Separate classical
worldsheet continuation from possible **quantum annihilation**.

Stage E: **physical recoil**. Derive endpoint momentum and charge
impulses in the conformal massive-boundary equations for pion,
electron+neutrino and muon+neutrino emissions; enforce four-momentum
and gauge constraint matching, test finite signal propagation
from the endpoint into the string, and rerun convergence.

## Decision gate

A full 1000 ps 3D classical worldsheet trajectory is **NOT VALIDATED**
until stages A–D are independently passed.
Radiation/annihilation assumptions and any collider reach remain
out of scope until stage E and the survival model are validated.

**Do not merge this research branch to main.**
