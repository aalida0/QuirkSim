# QuirkPheno AF archive audit — 2026-11-11

Verdict: execution smoke test passes; outbound zero-field analytic check
supported; 2 T accuracy under refinement remains unresolved. No recoil,
annihilation or solver replacement claim.

Input SHA-256: 5cba77d6e0c64bea5cd1fdd69054a7c736bd32a9b8ad87c9d7f1fddaca7b1432
Source: d91cb27ed3786bc63900faad892513841b5b6ed4.
Archive manifest and build log confirm Geant4 11.4.1 with LCG_110.
All macro hashes match. Build and four executions exit zero without timeout.
One event and both endpoints present in each ROOT file. Per-endpoint saved
step counts equal event-reported counts; times increase strictly and numerical
values are finite. No G4Exception marker found in case logs.

| Field | Level | End times, endpoint 0 / 1 (ps) | Step rows |
|---|---|---|---|
| 0 T | 0 | 21.4906295813 / 20.8686256297 | 17 / 18 |
| 0 T | 1 | 20.6827207647 / 20.8685578989 | 25 / 26 |
| 2 T | 0 | 21.4906295811 / 20.8686199208 | 17 / 18 |
| 2 T | 1 | 20.6827150544 / 20.8685578987 | 25 / 26 |

All eight terminal reason codes are 6, kTime in upstream Actions.hh, and
cover 20 ps. The stopping condition acts after exceeding a step, explaining
the overshoot. Final positions must not be compared at unequal final times.

## Independent analytic check

Use the ROOT-observed initial momentum, energy and mass, not the unrounded
macro numbers. In the zero-field outward branch of the symmetric straight
string, boost along z with b=pz0/E0 and g=1/sqrt(1-b^2). With tension F:

p_abs(t)=|px0|-F*c*t/g;
E_star(t)=sqrt(m^2+p_abs(t)^2);
x(t)=sign(px0)*(E0/g-E_star(t))/F;
z(t)=b*c*t;
pz(t)=g*b*E_star(t).

Use c=0.299792458 mm/ps and hbar*c=1.973269804e-10 MeV mm.
F=3035.955847424501 MeV/mm. Compare at each saved time, no interpolation.
Largest position residual over both B=0 runs/endpoints: 6.441514e-11 mm.
Largest momentum residual norm below 5.8e-6 MeV. These are selected
short-outbound analytic checks, not a general validation certificate.
The prediction is first turnaround at approximately 111.4298 ps;
these 20 ps runs do not exercise turnaround or string crossing.

## Refinement and magnetic signal

Compare over 0–20 ps using 2,001 common-time samples and cubic Hermite
position interpolation with endpoint velocities from ROOT momenta.
Confinement reference hbar*c/Lambda = 2.54944419e-7 mm.

| Field | Endpoint | Max coarse/fine position difference (mm) | Reference lengths |
|---|---|---|---|
| 0 T | 0 | 2.168412e-9 | 0.00850543 |
| 0 T | 1 | 2.322714e-9 | 0.00911067 |
| 2 T | 0 | 9.200962e-7 | 3.60901 |
| 2 T | 1 | 9.042126e-7 | 3.54670 |

At 20 ps, coarse/fine y deflections are -13.81584/-14.71585 nm for
endpoint 0 and -13.83274/-14.73695 nm for endpoint 1.
The approximately 0.9 nm resolution change is about 6% of the finer
deflection. It is small relative to the millimetre excursion but material
for effects on the confinement-length scale.

Interpolation has not been bounded independently for B=2; therefore these
are exploratory comparisons, not certified solver errors. Linear interpolation
gives approximately 8.6e-4 mm coarse/fine differences even in the analytically
accurate B=0 case. The B=0 analytic check demonstrates why that interpolation
must not be mistaken for a transport defect.

Event CPU is 0.00128–0.00150 s; this excludes initialization/build and is not
a speed comparison to the long, decaying conformal runs.

## Decision

Public-code reuse remains promising. Next establish magnetic refinement at
more levels and extend a bounded stable test through a turnaround. Do not
call the 2 T configuration validated for recoil yet. No repeat of the doublet
runs is needed. A shared initial string state and the charge-changing weak
decay adapter remain required for the intended recoil comparison.
The archive does not expose complete string energy-momentum histories; no
total-system conservation assertion can be made from these ROOT endpoints.

Reproduce with numpy, scipy and uproot:

    python scripts/audit_quirkpheno_transport.py quirkpheno_transport_pilot.tar.gz --output audit.json

Derived numbers and scope: studies/data/quirkpheno_transport_audit_2026-11-11.json.
