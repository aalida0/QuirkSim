# Triplet first-return window: selected classical recoil effect

The uploaded six-case AF pilot extends to 260 ps. The uploaded compact NPZ
archive contains steps between 220 and 226 ps. Its SHA-256 is stored with the
[reproducible numerical audit](data/triplet_first_return_window_2026-11-11.json).
The backend source is QuirkPheno d91cb27ed3786bc63900faad892513841b5b6ed4,
with the experimental localized weak-handoff adapter. All six macros match
their manifest hashes; all six processes exited zero with endpoint time-limit
codes. These are numerical transport observations, not annihilation events.

## Results

| Case | Level | Piecewise-linear minimum (micrometres) | Minimum-separation lower bound over the common uploaded window (micrometres) |
| --- | ---: | ---: | ---: |
| Stable | 7 | effectively zero | 0 |
| Stable | 8 | effectively zero | 0 |
| Sampled recoil, seed 1 | 7 | 4.154353 | 0.597009 |
| Sampled recoil, seed 1 | 8 | 4.154395 | 1.904739 |
| Opposite recoil, seed 1 | 7 | 4.028178 | 0.363428 |
| Opposite recoil, seed 1 | 8 | 4.028210 | 2.336632 |

The level-8 estimated minimum times are 223.015133 ps (stable),
222.872035 ps (sampled) and 222.994240 ps (opposite). The two recoil minimum
estimates shift by 0.0414 and 0.0313 nm between levels 7 and 8. Agreement
between these two levels is an empirical refinement comparison; it does not
bound errors in the underlying transport model.

The matching radius 8 hbar*c/Lambda is 2.039555 nm. The recoil cases remain
outside it throughout their common uploaded windows under the stated
kinematic assumptions, without interpolating to establish this exclusion.
The stable sample admits a minimum no larger than approximately 3.868 nm
using only the saved positions and speed-of-light bound; its data alone do
not establish entry into the 2.040 nm surface. Its interpolated zero must not
be interpreted as femtometre-scale numerical precision.

The earlier coarse monitor's first unresolved intervals were not observed
contacts. The stronger analysis resolves those ambiguities for both recoil
directions in this window.

## Bound and validation

For recorded endpoint positions x_a(t_a) and x_b(t_b), at time t the actual
separation is at least |x_a-x_b| - c(|t-t_a|+|t-t_b|), assuming continuous
subluminal endpoint histories. Take the largest bound from the four pairs of
adjacent observations. The actual separation is 2c-Lipschitz. With grid
spacing h, subtracting c*h from the smallest grid lower bound covers the
whole continuous time window. A 1e-9 mm engineering roundoff margin is also
subtracted. This is a mathematical kinematic bound with a floating-point
guard, not a formally certified interval-arithmetic result or a bound on the
Geant4 integration error. It remains conditional on the represented endpoint
histories, including continuity at handoffs. Only the overlapping span of
the two uploaded endpoint histories is covered, not arbitrary times outside
the archive.

Seven close-return tests include a known straight-line miss, an exact
crossing, asynchronous trajectories, invalid handoffs and 10000 randomized
bound checks. Inspection of the full-run logs gives eight preserved
handoffs, maximum four-vector component residual 5.66e-11 MeV, maximum
on-shell energy residual 5.83e-11 MeV, and maximum proper-clock root residual
8.34e-10 ps. The ROOT files were not included in the full-run archive; this
audit uses their separately supplied NPZ window records.

## Physical scope and release

This establishes a selected classical recoil-induced miss of the first
return in the specified forced-pion, fixed-kinematics seed. It is not an
ensemble probability or a physical quantum annihilation suppression factor.
Both endpoints have already undergone weak decay before the window; the
charged compact-state partial widths calculated separately are not applied
to these neutral daughter pairs. Quantum-state matching, neutral-sector
annihilation amplitudes and radiation/capture remain unresolved. No
annihilation event has been inferred or sampled from the closest approach.

Reproduce with:

```bash
python scripts/audit_quirkpheno_return_window.py triplet_return_window.tar.gz \
  --output triplet_first_return_window.json
```
