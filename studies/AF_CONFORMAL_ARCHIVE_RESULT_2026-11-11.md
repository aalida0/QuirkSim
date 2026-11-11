# Two-resolution AF trajectory audit

Input: user-supplied `conformal_convergence_diagnostics.tar.gz` from the eight
AF jobs at `e1ea0e08aa44acab15d68cadcd4ebb8aea12542a`.
SHA-256: `78313a66673301b20c1e638ddf0a54db0b5a9746838998a4eb6c63df82f1401e`.
No solver rerun or physical survival calculation was performed for this audit.

## Result

Both triplet seeds pass the existing numerical gates at both resolutions and
reach 1,000 ps. The saved trajectories also agree at less than 0.1 confinement
reference lengths. This is a two-resolution consistency result for these
initial conditions, not an established convergence order or validation of
quantum annihilation.

All doublet runs fail. They violate the closure tolerance before the eventual
stopping condition. More wall time or a larger step budget alone cannot
establish validity. No new production or tighter batch is warranted yet.

| Model, seed | Level 0 stop (ps) | Level 1 stop (ps) | Maximum endpoint difference at matched saved times / reference length |
|---|---:|---:|---:|
| Triplet, 1 | 1000, target | 1000, target | 0.0615442 |
| Triplet, 25 | 1000, target | 1000, target | 0.00603492 |
| Doublet, 1 | 315.081016, floating-point label | 262.077349, step budget | 0.928137 (shared interval only) |
| Doublet, 25 | 341.773919, step budget | 315.197447, step budget | 5.46423 (shared interval only) |

The reference is the submitted matching radius divided by eight: approximately
2.54944e-8 cm for the triplet and 9.16096e-9 cm for the doublet. Comparisons use
nearest saved observation times differing by at most 1e-7 ps; no interpolation.
The actual largest mismatch is 1.024e-11 ps, whose single-endpoint speed bound
is below 3.1e-13 cm, much smaller than the worst recorded position differences.
There are 218 matched records for each triplet seed, 61 for doublet seed 1,
and 78 for doublet seed 25. Dense near-decay records are not independent tests.

Triplet endpoint differences at exactly 1,000 ps are 1.56904e-9 cm (seed 1)
and 1.53857e-10 cm (seed 25), also the largest saved position differences.
Endpoint momentum differences there are 6.29447e-8 and 1.97001e-8 GeV.
Weak endpoint sequences agree, and maximum weak-decay time shifts are
1.012e-11 and 9.06e-14 ps. Reflection counts agree at 155 and 219 respectively;
this does not assert equality of every reflected ray or quantum phase.

## Doublet failure precedes termination

| Seed | First saved closure failure, level 0 (ps) | Level 1 (ps) |
|---|---:|---:|
| 1 | 78.235787 | 108.235787 |
| 25 | 162.149343 | 242.149343 |

These are first failures in the saved observations, not exact threshold-crossing
times. Fine runs improve the closure residual but terminate earlier than their
coarse counterparts; final residuals must not be compared as equal-time results.
Total energy conservation remains excellent (roughly 3e-11 to 8e-11 GeV
residual), demonstrating why energy alone is an insufficient numerical gate.

At the doublet stops, endpoint separations range from 2.16e-4 to 4.26e-3 cm,
far above the 7.33e-8 cm matching radius. The energy-coordinate intervals are
0.00820 to 0.14210 ps, far above the explicit 1e-8 ps chart cutoff. Neither
these endpoint values nor the stop labels establish a physical contact or
annihilation event. Saved five-ps observations cannot exclude unresolved
intermediate encounter minima independently of the solver's event detector.

## What the evidence does and does not identify

Source inspection shows that `max_steps=200000` is a budget for each call to
`advance`, not the total job. Calls normally cover at most five ps. The three
budget stops therefore reflect severe local stepping cost. The floating-point
label stop means a proposed outgoing null coordinate did not increase at
machine precision. It does not mean the physical string has annihilated.

The source combines RK4 endpoint evolution, spherically interpolated ray
histories and event-root handling at reflected weak discontinuities. The
saved files do not contain the last rejected trial, step size, rejection count,
active event caps or last null-coordinate increments. Consequently they do
not distinguish rapid physical ray rotation from repeated event truncation
or interpolation error as the cause of the stalls. Assigning one of these
as the proven root cause would exceed the evidence.

The next solver change should add those stop diagnostics and a targeted
short-interval reproducer, then improve accuracy control for ray histories
and event stepping as indicated by that evidence. Keep failed outputs and the
0.1-reference closure threshold. Do not project endpoints onto a reconstructed
chord, relabel numerical stops as quantum matching, or increase tolerances
to manufacture a pass.

## Reproduction and scope

```bash
python scripts/audit_conformal_archive.py conformal_convergence_diagnostics.tar.gz \
  --output studies/data/af_conformal_archive_audit_2026-11-11.json
```

The audit uses the standard library and reads archive members without extracting
them. The saved JSON contains provenance, stop times, first sampled closure
failures, paired comparisons and limitations. A separate NumPy calculation
independently reproduced the distance and timing summaries.

Verdict: selected triplet classical consistency supported; doublet long-time
numerical validation blocked. The radiation-recoil kernel, quantum matching
population and complete annihilation matrix remain physically unvalidated.
No unconditional annihilation-suppression or collider claim is enabled.
