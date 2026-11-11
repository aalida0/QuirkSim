# QuirkPheno 130 ps refinement result — 2026-11-11

Verdict: all eight execution smoke checks pass and the zero-field trajectory
agrees with the radial analytic solution through first turnaround. Magnetic
refinement improves systematically, but confinement-scale trajectory accuracy
is not yet established. No weak recoil or survival result.

Archive SHA-256: 92ec3b7af67ee64af93a28afde3daafb203c7a40df0e70b17d2dddde8b34151b
Pinned upstream: d91cb27ed3786bc63900faad892513841b5b6ed4.
Mass=225 GeV stable control, Lambda=774 eV, vacuum, beta=.4 in pair CM,
CM beta_z=.25, Bz=0/2 T. It is not the loop-split charged triplet parent
or the stretched initial worldsheet of the earlier conformal pilot.

## Execution and turnaround

Build and all eight runs exit zero without timeouts; all macro hashes match.
Each file contains one event, both endpoints, finite monotonically timed step
records, saved counts matching reported counts, and terminal reason 6=kTime.
All endpoint end times exceed 130 ps, ranging 130.10344–131.13510 ps.
No G4Exception marker found in case logs.

All sixteen endpoint trajectories change px sign. Their saved brackets include
the approximately 111.43 ps radial prediction; this is a bracketed sign change,
not an accurately solved magnetic turnaround root. At level 4 the brackets are:

| Field | Endpoint | First saved px-sign-change bracket (ps) |
|---|---|---|
| 0 T | 0 | 111.358768–111.472943 |
| 0 T | 1 | 111.351402–111.465610 |
| 2 T | 0 | 111.358773–111.472943 |
| 2 T | 1 | 111.351892–111.466098 |

Zero-field maximum analytic position residual is 2.842862e-9 mm,
approximately 0.01115 confinement reference lengths. Momentum residual norm
is at most 3.474e-5 MeV. These selected endpoint checks support transport through
the first turnaround; they do not measure total endpoint-plus-string invariants.
The analytic solution and constants are given in the earlier 20 ps archive
report and reproducible audit script.

## Magnetic refinement

Cubic Hermite position interpolation uses ROOT endpoint velocities, over 2001
common times from 0 to 130 ps. Worst endpoint difference for each adjacent pair:

| Levels | Max position difference (nm) | Confinement reference lengths |
|---|---|---|
| 1–2 | 1.40423 | 5.50797 |
| 2–3 | 0.705412 | 2.76692 |
| 3–4 | 0.350302 | 1.37403 |

Differences approximately halve on each refinement. This is evidence of
approximately first-order combined resolution sensitivity for these settings,
not a proven asymptotic error order. maxBoost and maxMerge were halved together,
so a unique error source cannot be assigned from this scan.

At 130 ps, the finest transverse y displacement is approximately -439.808 nm
for both endpoints. The level-3/4 displacement change is approximately 0.341 nm,
or 0.078% of that deflection. Small fractional change relative to the macroscopic
deflection does not establish accuracy for confinement-scale encounter physics.
Reference length hbar*c/Lambda=0.254944419 nm.

Transverse linear-vs-Hermite sensitivity (13001 common times) at level 3 is
at most 0.001474 nm and at level 4 at most 0.0003692 nm, much smaller than the
0.35 nm magnetic resolution difference. This supports interpreting the dominant
transverse difference as resolution sensitivity rather than merely sparse output
interpolation. It is not a rigorous interpolation error bound. Full-position
linear interpolation has much larger x errors and is unsuitable as the primary
comparison.

Event CPU per case ranges 0.00191–0.00717 s, excluding Geant4 initialization,
build and file staging. No speed comparison with decaying long-time runs
or projected production throughput is justified from that timing.

## What follows

The public backend is a credible reuse candidate for stable triplet-scale
transport. Four more halvings would put the next adjacent difference below
0.1 reference lengths if the observed trend persists; that is a forecast,
not a pass. A bounded finest-level confirmation can test it cheaply rather
than restarting custom doublet development.

The weak-decay/string-history adapter can be developed in parallel conceptually,
but no physical recoil result can be certified yet. Preserve upstream string
and partner state through the transition; do not activate crossing restart to
make a decay work. Actual crossing near the zero-field 222.86 ps re-encounter,
quantum matching, radiation and annihilation remain untested by these runs.
No further AF job was submitted by this audit.

Reproduction:
    python scripts/audit_quirkpheno_transport.py quirkpheno_refinement.tar.gz --output audit.json

Machine-readable result:
studies/data/quirkpheno_refinement_audit_2026-11-11.json
Transverse interpolation diagnostic:
studies/data/quirkpheno_refinement_interpolation_2026-11-11.json
