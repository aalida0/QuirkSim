# Frozen neutral-triplet phenomenological closure

## Deliverable and release

The radiation/action/capture/annihilation procedure is now explicitly defined,
implemented and evaluated on all 42 fine neutral-qualified fits. It supports
conditional phenomenological tables. It does not determine a unique physical
inclusive annihilation rate. The spin occupation, off-gate radiation and
matching probability are not inferred by naming a prescription.

## Fixed assumptions and justification

| Component | Adopted rule | Evidence / qualification |
|---|---|---|
| Model | Dirac Y=0 EW triplet, IC fundamental, N_IC=2, m0=225 GeV, Lambda_tension=774 eV, alpha_IC=.04 | Existing fixed benchmark, both-neutral vacuum evolution after first return |
| Energy radiation | Mean epsilon*Lambda per full signed radial reference period; epsilon=0,.01,.1,1 | Li et al. 2311.15486v2 Eq.(7) motivates mean loss; period convention is explicitly fixed, not silently transferred |
| Emission geometry | r_peri<=g/Lambda, g=1 or10 | Kang/Luty 0805.4642v3 Sec.5.7 close-approach argument; factor10 is a declared radius sensitivity, not a source error bar |
| Angular evolution | Conserved L or maximal downward one hbar per effective Lambda emission | Envelopes of bounded-step prescriptions; maximal drift is not the literature's unbiased randomization law |
| Capture | Absorbing boundary E<=8 Lambda and L<=1, matching efficiency=1 | Idealized matching prescription; selected thresholds E=1,32 Lambda and L=0,10 were stress-tested; overlap not calculated |
| Spin | No default singlet fraction | Vector-current production does not justify statistical 1/4; retain unknown sector explicitly |
| Singlet annihilation | Existing specified compact-state width 1.9750047986576727e-5 GeV | LO coupled-channel spin-zero calculation; lifetime 3.33e-20 s once populated |
| Other spin sector | Unknown annihilation probability between zero and one | No reuse of the singlet width for a triplet state |
| High-L hard annihilation | Neglected in the scenario tables, with a central WKB barrier diagnostic | Kang/Luty Eq.(5.27) motivates suppression; the diagnostic assumes a narrow high-L semiclassical event and is not a measured partial-wave distribution |
| Additional excited hazard | Zero in the frozen family | Previous .01 hazards remain stress tests, not derived rates |

Primary sources: https://arxiv.org/html/2311.15486v2 ;
https://arxiv.org/pdf/0805.4642 ; https://arxiv.org/html/2504.02940v2 .
The last source supports retaining deexcitation model dependence; scalar hard
annihilation widths are not transferred. The scale called Lambda here is
sqrt(tension); other papers' confinement conventions are not numerically equated.

The radial clock and absorbing transition are phenomenological reductions.
Linear-clock evolution is not asserted to resolve the final Coulomb cascade.
The cooling integral is dominated by the initial tens of GeV; threshold
variations are diagnostics of that approximation, not a determination of the
capture overlap. No sampled radiation four-momenta or NG backreaction are added.

## Correct spin treatment

For captured probability C, age t and specified singlet fraction f, the known
annihilation contribution is A=C*f*(1-exp(-Gamma0*t/hbar)). The uncalculated
sector gives the conditional interval [A,A+C*(1-f)]. If f is unspecified,
retain [0,C]. At exactly zero age all finite-width post-capture decays vanish.

Therefore the previous f=1/4 result of 25% is the known singlet contribution
at late times, not a prediction that total annihilation is exactly 25%.
Its interval is [25%,100%] if capture is certain. Full occupation f=1 remains
a prepared-singlet scenario, not the DY production default.

## Benchmark findings

* All 14 isotropic recoil fits and 22 collinear 2 T fits have zero capture
  in the tested frozen family. Their reported zero annihilation is conditional
  on this no-excited-hazard, contact-gated prescription, not an absolute zero.
* The six zero-field collinear controls capture in the efficient epsilon=.1
  baseline at about .151 seconds. Their total annihilation interval is [0,1]
  without a specified spin population; it becomes [1,1] at late times for the
  explicitly prepared singlet. Epsilon=0 remains a no-cooling alternative.
* No tested fine case captures within one nanosecond of the fitted handoff.
  Thus these results do not establish a detector-time annihilation contrast.
* This is a both-neutral subset, not an inclusive event survival estimate.

## High-L cross-check

For the central Hamiltonian H=2sqrt(m^2+p_r^2+L^2/r^2)+sigma*r-kappa/r,
take 0<=kappa<=.2. This ceiling exceeds the fixed benchmark short-distance
attraction CF*alpha_IC+2*alpha_W, about .098. In the checked domain
1/m<=r<=L/H, L>=2*kappa and sigma*L<H^2, the forbidden radial momentum
satisfies |p_r|>=sqrt(7)*L/(4r). Hence

    S >= sqrt(7)*L/4 * log(m*L/H).

We evaluate this after granting the maximum unit-step angular loss allowed
by the full kinetic-plus-binding energy budget, and using the initial energy
to reduce the barrier conservatively. For the isotropic fits, even the least
suppressed WKB exponential has log10(exp(-2S)) below -5.36e12. This explains
the negligible direct hard-annihilation approximation within the central
semiclassical model; it is not a rigorous physical probability bound including
unmeasured low-L tails, spin/flavor occupation or field modes.

## Validation and reproducibility

Sixteen focused tests pass. New checks preserve the unknown spin interval,
test zero capture/age and compare the analytic action floor with independent
radial quadrature. Saved-data assertions check all 42 short-horizon outcomes
and the high-L/zero-field-control interval classifications. Conservation
targets from the preceding audit support the fitted angular scale; they
do not replace direct string quadrature or quantum matching.

```bash
python scripts/finalize_triplet_closure.py --output frozen_closure.json
```

The saved output records source/stress/archive hashes, fixed prescription,
per-case spin intervals, emission-budget diagnostics, WKB action diagnostics
and release flags. It consumes existing validated numerical outputs and does
not run AF transport. Broad-gate/maximal-angular-loss counterexamples remain
in the separate stress artifact and are not removed from the evidence.

The implementation task is complete for this declared phenomenological
family. Claiming absolute physical closure would still require a justified
capture occupation/overlap and the uncalculated spin-sector rates, plus any
material off-gate or coherent radiation mechanism. Those are not supplied by
the inspected literature and have not been manufactured here.
