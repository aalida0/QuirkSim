# Radiation and capture prescription stress test

## Finding

The first-return geometry effect survives independently of capture modeling.
The annihilation contrast does not survive arbitrary radiation/action-loss
closures. The null additional suppression at 2 T is stable under a tenfold
larger emission gate, but is not a prescription-independent physical result.

Using the saved 42 level-8 neutral-qualified direction cases, 2352 secular
paths and 9408 paired closure cards were evaluated at a ten-second horizon.
No AF transport was run. The same compact-state width and binding inputs as
the previous comparison were retained.

## Interventions

* Emission radius: 1, 10, 1000 and 100000 times 1/Lambda.
* Mean energy emission: epsilon = 0, .01, .1, 1 per full reference period.
* Directed angular loss: 0, 1, or the gate multiple per effective Lambda emission.
* Capture threshold: default 8 Lambda and L<=1, plus selected cards at
  1 or 32 Lambda and L<=0 or L<=10 (not a full Cartesian uncertainty scan).
* Capture singlet occupation: 0, 1/4, 1; excited optical depth per eligible
  encounter: 0, .01. These are separate competing-outcome assumptions.

Large gates and maximal directed loss are stress interventions, not published
emission kernels or calibrated nonperturbative uncertainties. In particular,
the angular bound built into the prescription permits loss up to the gate
multiple, but does not derive such directed emission dynamically.

| Gate multiple | 2 T paired isotropic-minus-control probability range over valid tested cards |
|---:|---:|
| 1 | 0 |
| 10 | 0 |
| 1000 | -1 to 0 |
| 100000 | 0 |

These are individual paired-card extrema, not mean contrasts or confidence
intervals. The final row does not imply universal agreement: 240 paths at
that gate leave the allowed central-orbit domain, so their outcomes remain
unknown. Valid broad-gate, maximal-loss cards can capture both direction
modes, erasing the contrast. Intermediate gates can admit collinear orbits
while leaving isotropic recoil outside, creating a negative contrast.

## Independent budget check

For energy loss in units of Lambda and constant directed loss coefficient a,
a necessary capture condition is

    a >= max(0,L_initial-L_capture)*Lambda/(E_initial-E_capture).

This follows from the integrated angular and energy ledgers and is necessary,
not sufficient: the gate can close or the central reduction can fail first.
For isotropic cases the required coefficient spans about 6823--51165.
The collinear controls at 2 T also require coefficients far above one;
zero-field collinear controls have symmetry-protected L=0. Thus simply
opening the radiation gate cannot make the original one-unit-loss model
capture these large-L recoil states within its energy budget.

Initial isotropic gates open only at radius multiples about 2849--25925.
Collinear controls span zero at 0 T to about 1564 at 2 T. This explains why
factor-ten gate variations leave the original classification unchanged.

## Validation and claim decision

Energy/angular ledgers were checked across every path. Four representative
benchmark paths were repeated with half the energy step; their status and
evolution times agree to the saved tolerances. Existing focused evolution
tests were also run. These are internal/model checks, not external physical
validation. Details are in `data/triplet_capture_stress_checks_2026-11-12.json`.

The full cards, including truncated outcomes, are losslessly saved in
`data/triplet_capture_stress_2026-11-12.json.gz`. The input artifact hash and
original AF archive hash are recorded. The script reproduces the cards:

```bash
python scripts/stress_triplet_capture.py \
  studies/data/triplet_ensemble_capture_models_2026-11-11.json.gz \
  --output stress.json.gz
```

Release decision: report the geometric recoil effect and conditional
annihilation scenarios. Do not report a unique physical annihilation
probability or use the card range as an uncertainty band. A credible
annihilation prediction needs a justified relation between energy radiation
and angular-momentum transport, including any off-gate emission; more seeds
alone cannot determine that relation. The present excited hazard is still
gated, not an implementation of off-gate radiation. Mixed-charge encounters,
pre-handoff losses, string modes and quantum spin occupation remain outside
this neutral secular reduction.
