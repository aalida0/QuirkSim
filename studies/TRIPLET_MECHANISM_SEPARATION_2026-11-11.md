# First-return mechanism comparison

Step 1 is complete using the existing 116-case compact AF archive. No new
transport was run. The script independently recomputes all CM chord minima
from NPZ windows; all agree with the archived broad-search minima. It saves
all decay categories, paired seed differences, field contrasts and refinement
checks in `data/triplet_mechanism_comparison_2026-11-11.json`.

## Results

The table uses level 8 and the same both-neutral seed subset in each of the
three decay-angle modes. Stable transport is one deterministic control per
point, shared by the seed comparisons, not an additional random sample.

| Point | Qualified seeds | Stable fitted miss (nm) | Collinear forward/backward median miss (nm) | Isotropic median miss (nm) | Median paired isotropic-minus-forward miss (nm) |
|---|---:|---:|---:|---:|---:|
| 2 T, beta .4, theta 90 | 3 | <0.001 | 70.67 / 71.62 | 4188.43 | 4117.76 |
| 2 T, beta .4, theta 60 | 3 | 0.0201 | 61.20 / 62.02 | 3625.13 | 3563.93 |
| 2 T, beta .6, theta 90 | 5 | <0.001 | 188.33 / 189.50 | 4870.24 | 4668.41 |
| 0 T, beta .4, theta 90 | 3 | <0.001 | consistent with zero | 4127.32 | 4127.32 |

Differences are computed per matched seed before summarizing. Subtracting
separate medians is generally different. Distances are nonlinear scalar
observables, so these contrasts are not an additive force decomposition.

Magnetic bending of the stable charged pair does not by itself produce a
large fitted first-return miss at these points. The collinear charge-changing
decays introduce a finite miss at 2 T, while their zero-field counterparts
remain consistent with a one-dimensional crossing. Allowing transverse pion
recoil then increases the miss to micrometre scales in both fields.

Stable-to-collinear changes mass, longitudinal momentum and charge history
together. Consequently this comparison identifies a combined decay/magnetic
baseline, not a pure charge-history effect. Collinear-to-isotropic compares
conditional decay-angle interventions with common proper thresholds. The
benchmark/zero-field cases share initial kinematics and latent plans, but
field-dependent evolution changes the actual decay vertices.

## Checks and limits

All nine zero-decay level-8 cases reproduce the corresponding stable minimum
exactly. One-decay categories also exhibit the hierarchy: isotropic median
misses are 3.46--5.66 micrometres; collinear medians at 2 T are 93--311 nm,
while zero-field collinear values are consistent with zero. These are geometry
checks only; no mixed charged-neutral capture prescription is applied.

The largest stable level-7/8 shift is 1.43e-5 nm. The previously recorded
full mode refinement maximum is 0.047 nm. Interpolated minima do not certify
contact: speed-only bounds remain unresolved for the stable 2 T cases.
Zero-field stable contact within the numerical matching radius is independently
certified by those bounds. The 0.255 nm emission gate and 2.040 nm matching
radius must not be conflated.

This resolves the mechanism-comparison task, not inclusive annihilation.
The existing central capture model still gives no extra annihilation contrast
between collinear and isotropic samples at these 2 T points, despite their
large geometric difference. Radiation/gate robustness is the next question.
Eight shared seeds, forced pion decays and no production weights limit the
population interpretation.

Source archive SHA256:
`ca6ff21c5447efa156a65732430be467d27ba993ab8218558e329d6be6773cf1`.

Reproduce:

```bash
python scripts/separate_triplet_mechanisms.py triplet_ensemble_v1_audit.tar.gz \
  --output mechanism_comparison.json
```
