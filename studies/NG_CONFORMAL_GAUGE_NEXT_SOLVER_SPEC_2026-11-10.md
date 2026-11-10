# Step 1 next solver: conformal/double-null Nambu–Goto with massive charged boundaries

**Date:** 2026-11-10  
**Status:** PARTIALLY IMPLEMENTED: zero-field boundary flux and causal
characteristic prototype, analytic radial/rotating checks and pulse refinement.
See [implementation and remaining gates](CONFORMAL_BOUNDARY_RESULT_2026-11-10.md).
No validated general contact continuation or physical recoil is established.

## Why change coordinates?

The zero-artificial-mass equal-spatial-arclength gauge reproduces the
exact radial benchmark and converges pre-contact at N=15–25, but stalls
when a boosted 2T string contracts to a few confinement lengths.
Crucially, last-valid diagnostic at N=11 and total arclength
2.823751e-7 cm gives:

- fastest **vertex-label** speed `0.999993799 c`;
- fastest **transverse midpoint string segment** speed `0.932873923 c`;
- minimum segment Nambu–Goto lapse-squared
  `1-v_perp^2 = 0.129746243`.

Hence the failure is NOT proof of a physical lightlike cusp. It is
a strong reason to stop using a spatially equidistributed gauge
through severe contraction and folding.

At requested total length 1e-7 cm, DOP853 and Radau both encounter
the minimum floating-point integration step around the same length
(~2.82e-7 cm). BDF fails earlier (~5.49e-7 cm) and LSODA reports
a numerical superluminal vertex. These are **fail-closed** outcomes.

## Formulation to derive before writing the full integrator

Use a covariant classical NG action with massive dynamical endpoints
and external electromagnetic boundary coupling:

```text
S = -T ∫_Σ d²ξ sqrt(-det h_ab)
    -Σ_i m_i ∫ ds_i
    +Σ_i q_i ∫_{∂Σ_i} A_μ dX_i^μ,

h_ab = ∂_a X^μ ∂_b X_μ,
T = Lambda_IC²   (in natural units).
```

In flat spacetime, a conformal gauge gives the bulk linear
worldsheet equation (natural units):

```text
∂_τ² X^μ - ∂_σ² X^μ = 0,
(∂_τ X ± ∂_σ X)² = 0.
```

If additionally an orthonormal/temporal gauge is valid
(`X^0=t` and a suitable energy coordinate `σ`),
the spatial constraints become

```text
|∂_τ x|² + |∂_σ x|² = 1,
∂_τ x · ∂_σ x = 0.
```

**Warning:** Massive endpoint boundary coordinates generally must
be treated as moving `σ_i(τ)` unless a consistent residual
reparameterization is established. Do NOT freeze the endpoint
labels and impose massless Neumann conditions. Boundary
variation must reproduce both the massive endpoint momentum
equation and the string energy–momentum flux, including the
electromagnetic Lorentz force.

A suitable covariant boundary form is conservation of endpoint
four-momentum against outgoing NG canonical worldsheet flux
and `q_i F^{μν} u_ν`; the detailed normal orientation, moving
boundary terms and proper-time factors must be derived from
the action rather than guessed.

## Independent source checks

1. [Goustidis, Pourtsidou, Sakellariadou, Phys. Rev. D 91,
   025022 (2015)](https://doi.org/10.1103/PhysRevD.91.025022):
   conformal/temporal string gauge, moving massive junction
   and proper variation of boundary conditions. A massive
   endpoint is a one-leg specialization and must be rederived.
2. [Global conformal gauge in string theory,
   Phys. Lett. B 2021](https://doi.org/10.1016/j.physletb.2021.136246):
   treatment of conformal gauge and degeneracy of induced worldsheet
   metric at open-string boundaries.
3. [Sonnenschein, Weissman, Yankielowicz, JHEP07(2020)156](https://doi.org/10.1007/JHEP07(2020)156):
   classical rotating folded open strings with opposite
   electric charges in a constant field; mass and
   phenomenological setup differ from the quirk benchmark.
4. [Li et al., EPJC 80 (2020) 651](https://doi.org/10.1140/epjc/s10052-020-8209-y):
   approximate quirk endpoint equation including the NG
   transverse-velocity correction, magnetic field and ionization.
   Use as a *local boundary-force* benchmark, not a worldsheet solver.

## Implementation gates, in order

1. Independently derive and unit-test the flat-space conformal
   bulk and **massive endpoint boundary flux** equations in
   a consistent worldsheet coordinate system. Record all
   units and the sign of each endpoint current.
2. Prove that the selected coordinate/gauge constraints stay
   satisfied under a traveling transverse pulse and its
   reflection from a massive endpoint. Check causal propagation
   along the string and the time resolution of the pulse.
3. Reproduce the exact radial 1+1D massive-endpoint reference
   for B=0, including multiple signed crossings **only as
   a chosen analytic continuation**, not annihilation.
4. Reproduce the B=2T pre-pinch endpoint trajectories and
   energy accounting of the finite-element code before its
   failure; compare at the SAME physical laboratory time.
5. Perform a mesh and timestep refinement scan of **physical**
   observables, including endpoint chord and string energy
   at arclength ~4/Lambda (and smaller only when justified).
   Demand spatial uncertainty <0.1 hbar*c/Lambda before
   calling any short-distance geometry robust.
6. Do NOT infer annihilation at a geometric threshold.
   Only after a validated continuation or an explicitly
   parameterized absorbing boundary is derived should
   the model incorporate coherent quantum annihilation,
   radiative momentum backreaction, and weak-decay impulses.
7. Preserve early-contact survival: annihilation/radiation
   before weak decay must not be ignored in conditional samples.

## Decision gate for the pheno paper

**Not yet GO.** A compliant conformal/double-null solver could fix
the current gauge problem, but neither 1000 ps full 3D trajectories
nor an annihilation probability has been validated. No
LHC discovery/exclusion claim may use the present short-distance
geometric encounter fractions.

