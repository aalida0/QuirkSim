# AF triplet ensemble: direction controls change the annihilation interpretation

## Main finding

The 116-case AF pilot completed and passed its full on-AF local audit. The
uploaded archive is only 2,275,228 bytes because ROOT files stayed on AF.
Independent inspection of the packed logs confirms 174 four-momentum-
conserving decay vertices and preserved handoffs, including 116 exact
production-CM collinear control decays.

The transverse pion kick substantially increases first-return miss distances.
However, at 2 T the conditional collinear controls also miss the declared
capture gate. Thus the connected central model predicts **no additional
annihilation suppression from isotropic transverse recoil relative to these
controls at the three tested 2 T points**. The distinction is visible at zero
field: collinear controls approach the matching region and isotropic recoils
miss it. This is a narrower conclusion than one inferred from an artificial
zero-angular-momentum reference.

The collinear controls retain the mass changes, longitudinal recoil and
charge-history effects. They suppress the transverse kick in the production
CM and conserve local four-momentum. They are conditional decay-angle samples,
not an inclusive physical zero-recoil process.

## Geometry of the neutral-qualified level-8 subset

Each point has eight isotropic lifetime/direction seeds. The same seeds are
shared across points and controls. Only events with two neutral endpoints by
the estimated first return enter the current neutral capture model.

| Point | 0/1/2 decays by return, isotropic | Neutral-qualified seeds | Median CM miss, isotropic (micrometres) | Median CM miss, forward/backward collinear control (micrometres) |
| --- | --- | --- | ---: | ---: |
| beta=0.4, theta=90 degrees, 2 T | 1/4/3 | 1,4,7 | 4.1884 | 0.07067 / 0.07162 |
| beta=0.4, theta=60 degrees, 2 T | 1/4/3 | 1,4,7 | 3.6251 | 0.06120 / 0.06202 |
| beta=0.6, theta=90 degrees, 2 T | 0/3/5 | 1,3,4,7,8 | 4.8702 | 0.18833 / 0.18950 |
| beta=0.4, theta=90 degrees, 0 T | 1/4/3 | 1,4,7 | 4.1273 | consistent with zero / consistent with zero |

The 2 T collinear-control minima span roughly 42--399 nm, versus the
1/Lambda_tension gate of 0.255 nm and numerical matching radius 8/Lambda
of 2.040 nm. Seed-1 level-7/8 CM minima differ by at most 0.047 nm across
the tested points/modes. This is an empirical refinement check, not a formal
trajectory-error bound for every seed.

The speed-only continuous-time bounds independently exclude the 2.040 nm
radius in 9 of the 14 neutral-qualified isotropic point/seed cases. The
remaining five are unresolved by these conservative bounds. All zero-field
neutral-qualified collinear windows have an upper bound below about 0.001003
nm, including the engineering roundoff guard. At 2 T the collinear speed-only
bounds remain unresolved; their larger misses are interpolation/central-fit
results with the stated refinement evidence. An unresolved bound is not a
contact or a numerical transport failure.

## Numerical angular floor: corrected only by an exact control symmetry

Near a zero-field collinear crossing, tiny transverse numerical residuals
become hundreds of hbar after conversion to orbital angular momentum. Treating
this raw endpoint fit as physical would spuriously prevent the L<=1 capture
criterion even in a symmetry-protected one-dimensional control.

The analysis now sets L=0 **only** for the vacuum, zero-field, CM-collinear
controls after checking their transverse position and momentum residuals.
The raw L values and residuals are retained in the result. All eight such
windows (including refinement copies) pass, with maximum transverse position
residual 1.334e-12 mm and momentum residual 4.465e-11 GeV. Other modes and
nonzero-field data are not rounded or projected. A failed symmetry check
stops the fit. A new unit test checks both the permitted correction and
rejection of a genuinely noncollinear control.

The invariant follows from collinear initial conditions, zero external field,
a straight string and exactly CM-collinear decays. It is a calibration of a
declared symmetry limit, not a general procedure for assigning quantum L
from classical trajectories.

## Declared capture-model comparison

The 54 qualified windows include physical direction controls and seed-1
refinement copies. They generate 1296 secular orbit scenarios and 576 paired
conditional summary rows. Refinement copies are excluded from population
means. The current capture model is the explicitly declared central reduction
of `REPEATED_ENCOUNTER_RESULT_2026-11-11.md`.

For epsilon=0.1, directed angular loss, a capture singlet fraction 1/4,
zero excited-state annihilation hazard and a ten-second post-handoff horizon:

| Point | Mean isotropic-minus-collinear annihilation probability in the declared model |
| --- | ---: |
| Benchmark, 2 T | 0 |
| Tilted, 2 T | 0 |
| Faster, 2 T | 0 |
| Zero field | -0.25 |

Both forward and backward controls give these contrasts. Across the full
tested capture-spin/excited-hazard closures, the paired mean contrast is zero
at the three 2 T points and ranges from -1 to 0 at zero field. This is a
scenario range, not a statistical confidence interval or physical bound.
With a spin-preserving capture closure and no excited-state annihilation,
the contrast can vanish even at zero field. The zero paired standard errors
of saturated outcomes in some cards do not establish precision: there are
only three or five qualified seeds per point.

All 0/1-decay categories remain in the manifest and analysis summaries. This
subset result must not be converted into inclusive survival by discarding
those events. The eight common seeds across the four points do not form
32 statistically independent lifetime/direction samples. Production PDF
weights, inclusive non-pion decays, pre-handoff annihilation/radiation and
capture of mixed charged-neutral pairs are absent.

## Provenance and validation

The AF assets and the three recorded driver hashes match the checked local
source. The upstream source is d91cb27ed3786bc63900faad892513841b5b6ed4;
the new controls and wrapper were introduced by QuirkSim a40b333. All macro
and plan hashes match. The Geant4 build returned zero and the upstream
preimage checks passed.

Independent logged-vertex checks give maximum four-vector component residual
9.161e-11 MeV, mass-shell energy residual 5.821e-11 MeV and proper-threshold
residual 6.121e-10 ps. The maximum collinear rest-direction residual is
1.295e-15 in unit-vector norm. The full ROOT audit was performed on AF;
the ROOT files themselves were not uploaded and are not claimed re-audited
from this archive. Local checks here use their compact NPZ records and logs.

The source archive SHA-256 is
`ca6ff21c5447efa156a65732430be467d27ba993ab8218558e329d6be6773cf1`.
Saved outputs:

* `data/triplet_ensemble_summary_2026-11-11.json`: compact audit and comparison.
* `data/triplet_ensemble_window_bounds_2026-11-11.json.gz`: 116 window checks.
* `data/triplet_ensemble_capture_models_2026-11-11.json.gz`: full secular cards,
  angular diagnostics and conditional paired comparisons, compressed losslessly.

## Reviewer decision

The project now has an audited, physical conditional-angle control. It
supports a substantial transverse-recoil geometry effect and a zero-field
model contrast. It does **not** support a generic claim that transverse
recoil adds annihilation suppression in a 2 T detector field: the tested
collinear decay controls already block the gate in the current model.

The next scientific question is how the charge-changing/magnetic baseline
and the transverse kick act across the relevant parameter space. A paper
should separate those mechanisms and report model dependence. These small,
forced-pion, unweighted pilots and post-handoff capture closures are not yet
a production-level physical annihilation or detector-sensitivity result.

Reproduce from the uploaded compact archive:

```bash
python scripts/audit_quirkpheno_return_window.py triplet_ensemble_v1_audit.tar.gz \
  --output ensemble_bounds.json
python scripts/compare_repeated_encounters.py triplet_ensemble_v1_audit.tar.gz \
  --output ensemble_models.json
python scripts/summarize_triplet_ensemble.py triplet_ensemble_v1_audit.tar.gz \
  --bounds ensemble_bounds.json --models ensemble_models.json \
  --output ensemble_summary.json
```
