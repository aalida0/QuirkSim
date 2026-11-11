# Experimental public-backend weak-recoil pilot

This adapter is a patch against tvami/QuirkPheno commit
`d91cb27ed3786bc63900faad892513841b5b6ed4`. Upstream is Apache-2.0;
its NOTICE credits the ATLAS quirk extension and Geant4 transportation.
Keep upstream LICENSE/NOTICE in the checked-out source. This directory stores
an integration patch and new source, not a second independent transport solver.

**Status: first deterministic AF pilot compiled and local handoff checks passed.
The UI precision and localized clock passed the ten-case AF check at 66eef0b1.
Micrometre directional recoil is resolved in the deterministic pilot.
The new sampled-plan/two-endpoint extension is pending Geant4 verification.**
Do not merge into stable code.

The patch registers neutral Quirk definitions with sign-preserving endpoint
slots, transfers parent ownership through a one-use token, suppresses ordinary
endpoint death during that transfer, and bypasses production reinitialization
on the neutral daughter's first start. Fingerprints check both incoming string
histories, initialization/borrowed state and liveness across the handoff.
Cached endpoint velocity and kinetic energy are updated to the daughter.
No string restart, reset, projection or partner kick is introduced.

A post-step process performs a physical charged-parent -> neutral-daughter +
pion decay at proper time 50 ps. The default localized transport mode integrates
Gauss3 m/E dt on trial replay states and caps the accepted step at the event.
Legacy mode retains the first-accepted-step threshold rule as a control. This is a
**fixed-threshold pion integration test**, not a stochastic lifetime ensemble.
The computed clock root and the physical trajectory still need resolution review;
this is not an exact proper-time certificate for a continuum string solution. The parent neutral mass is 225 GeV; the parent charged
mass is 225.1571702763212 GeV from the existing QuirkSim LO triplet calculation.
The pion mass used in Geant4 is read from its particle table and logged.

The two rest-frame pion directions +/-y give physical angular controls. Neither
is a no-recoil limit. A stable-versus-decaying comparison mixes charge, mass,
and recoil. This pilot does not implement leptonic channels or all multiplet
production modes and does not certify the Dirac-triplet model globally.

### Close-return post-run censor

`scripts/audit_quirkpheno_close_return.py` checks every saved endpoint step in
single-event pilot ROOT files. It combines conservative speed-of-light bounds
over overlapping, asynchronous steps, including records on either side of a
weak handoff. It refuses missing steps and detects discontinuous sampled
handoffs. The monitor arms only after a separated interval, so its statement is
conditional on that arming; the point-production history before arming is not
covered. An unresolved interval is censored at its start; it is not evidence
of a contact, quantum transition, or annihilation. The rounding margin is an
engineering guard, not interval-arithmetic certification. As with all saved
step analyses, the speed and continuity assumptions apply to the represented
worldlines, not only their observed endpoints.

The pilot's optional `--audit-close-return` runs this audit after each case and
stores the compact result in the manifest. `--duration-ps` and
`--case-timeout-seconds` configure longer runs. The audit requires `uproot` in
the Python environment and `stepPrescale=1`. It is post-run: it does not stop
Geant4 online and does not restore a state at the censoring boundary. Full
trajectory and weak-handoff validation are still required for any physical
claim. Local checks: five synthetic tests and the existing 20 sampled-pion
pilot ROOT files. All 20 remained outside the matching radius over their
common histories from 5 ps through at least 130 ps, with weakest bound
1.127 mm; this is before the first-return region and is not survival through it.

## AF runner

`python3 scripts/af_quirkpheno_weak_pilot.py --output-dir <new-directory>`

Uses an isolated fresh upstream checkout under the output directory and applies
this patch plus the added headers/source. The earlier stable checkout is not
modified. Manifest records upstream revision, adapter file hashes, macros and
execution results. LCG_110, Geant4 11.4.1, one CPU, B=2 T, Lambda=774 eV,
point production, 130 ps. Eight localized cases at levels 7/8:
- stable charged parents;
- endpoint-0 pion +y and -y;
- endpoint-1 pion +y (charge-conjugate endpoint transfer).

Two additional level-8 controls use the same corrected inputs with the legacy
clock (stable parents and endpoint-0 +y pion), for ten cases in total.

No Condor. Per-case timeout 120 s (+5 s kill grace), build cap 600 s.
Runner stops on its first failed/incomplete case and keeps diagnostics.
All weak environment settings are explicitly assigned per case.
No sampled weak rate, radiation backreaction, annihilation, or survival claim.

## Evidence and local checks

`python3 research/quirkpheno_recoil/check_native.py`

Compiles native C++ with warnings as errors; checks local four-momentum,
on-shell daughter/products, closed channel rejection, one-use handoff and wrong
parent/daughter/history rejection, both endpoint signs. Eight canonical boost
fixtures agree with the independent existing Python pion sampler within
1.71e-13 GeV. It also checks the frozen charged mass against physics.benchmark.

With numpy/scipy/uproot/awkward:
`python3 research/quirkpheno_recoil/check_audit.py`

Builds synthetic ROOT/log fixtures, exercises the eight-case auditor including
piecewise interpolation across recoil, and verifies rejection of a corrupted
history-preserved record. Synthetic trajectories are not Geant4 evidence.
`git apply --check` succeeds against the downloaded pinned upstream snapshot.
The Python runners/auditor compile. No full Geant4 build was possible locally.

## Returned-data audit

`python3 scripts/audit_quirkpheno_weak_pilot.py <archive.tar.gz> --output audit.json`

Checks one decay/transfer in each weak case, local four-momentum and mass shells,
rest-frame pion direction, endpoint sign, neutral mass/charge, one event, track
coverage, macro identity, exceptions, row counts and history fingerprint record.
Records proper-time overshoot and physical-angle comparisons. Interpolates
separately before/after the recoil discontinuity. Parent terminal handoff rows
may carry zero kinetic energy after Geant4 kills the parent; logged pre-decay
four-vectors are the physical event state used instead of those rows.

Partner differences outside a conservative Euclidean future light cone are a
causality diagnostic, not a certified test of arrival along the actual string.
Full causal propagation, reflected fronts, total string energy-momentum and
refinement still require returned trajectories and additional checks.

## Patch preparation correction (2026-11-11)

The first AF invocation stopped at patch_check.log before compilation. The
original local snapshot had an extra final blank line; the appended StringStore
hunk therefore expected a line absent from the real checkout. This was reproduced
against byte-exact upstream files and corrected by regenerating the patch.
All eight preimage Git blob IDs now independently match GitHub's pinned tree.
The runner checks their SHA-256 values before applying and prints the end of
any failed preparation/build log. No whitespace-ignoring or fuzzy patch option
is used. Native and synthetic checks were rerun; full Geant4 build is still pending.

The next AF build exposed a scope error in the application instrumentation:
three added expressions in quirkStep used the caller-local name `def`.
They now obtain the definition from the local track `t`. The corrected patch
applies to the exact upstream bytes, and all three references were inspected
in the resulting function. Full executable compilation remains pending rerun.

## First returned weak run and UI precision correction

The d846f752 adapter compiled and all eight AF cases completed. The audit
reports one conserved pion event and one consumed handoff in each weak case.
Opposite physical directions change the decaying endpoint by 12.062 micrometres
and its partner by 0.700 micrometres at 130 ps. This is a direction comparison,
not a no-recoil or survival result; see the saved study for precision limits.

Actual backend mass logging exposed UI rounding of 225.157170... to 225.157 GeV.
The patch now enables G4UImanager::UseDoublePrecisionStr(true) before macro
commands, preserving input precision during Geant4 unit normalization. The
new main-file preimage is pinned and checked. That UI change and the localized clock subsequently passed the AF checks.
This establishes the computed-clock pilot, not continuum or quantum validity.

## Localized clock update

ProperTimeClock.h integrates m/E over lab time with three Gauss nodes, located
on replayed trajectory states before any string consumption. It brackets the
proper threshold in the curve-length parameter and caps the transport step.
The sampled replay cursor is restored to the accepted state by Update; the
sole consumption/push operation remains Update. Forced sub-navigation steps
now respect the event cap too. Localized mode requires moveTinySteps=true.

The computed root tolerance is 1e-9 ps, with 1e-12 ps lab-time node tolerance.
The event process accepts that numerical tolerance; the audit requires a
localized event residual within 2e-9 ps. These are tolerances on the computed
quadrature clock, not proven physical-clock errors. Legacy checks are unchanged.

Native C++ tests compare clocks and roots with exact constant-gamma and
constant-force references, including turnaround, and reject invalid brackets.
The synthetic audit covers all ten cases and the legacy-control comparisons.
The previously uploaded archive remains readable with unchanged local-check
results. The clock update subsequently compiled and passed the AF pilot; see the localized-clock result study.

## Sampled-lifetime forced-pion pilot

`python3 scripts/af_quirkpheno_weak_pilot.py --sampled-pion-seeds 8 --output-dir <new-directory>`

Prepares the complete seed block 1-8: independent exponential proper thresholds
for both endpoints and isotropic pion directions, paired with antipodal directions
at the same thresholds. Uses the existing nominal free-particle triplet lifetime
258.5200888747253 ps and reference pion branching 0.976048. The channel is forced
to pion: this is not an inclusive decay ensemble or a survival calculation.
The LO rate reference uses pion mass 139.57039 MeV; transport uses Geant4's
139.5701 MeV, and this small mass-scheme difference is explicitly retained.

Main seeds run at level 8; seed 1 also at level 7; two stable references give
20 cases total. Initial kinematics remain fixed, and the nominal window is
130 ps. Each process executes exactly one event; a worker's immutable plan
is not a multi-event schedule. `QP_WEAK_ENDPOINT=-1` enables both charged
endpoints; neutral daughters never receive another weak process.

Plans use a versioned, ordered two-row text format with times in ps and unit
rest-frame directions. Manifest records each file hash. The parser rejects
missing/reordered endpoints, nonpositive/nonfinite lifetimes, nonunit directions
and trailing data. Native tests and exact Python/C++ roundtrips pass.
`check_sampled.py` tests synthetic zero/one/two decays and bad handoffs; it
provides no Geant4 transport evidence. The new extension has not run on AF.

Audit with `scripts/audit_quirkpheno_sampled_pion.py`. Local checks include
proper thresholds, both endpoint handoffs, mass/charge/identity, physical pion
directions, four-momentum and censoring. Refinement/causal trajectory analysis
remains a separate returned-data task. Decay-count histograms refer to events
inside nominal 130 ps, not post-step end-time overshoots. Do not interpret them
as inclusive rates, collider efficiencies or quantum survival probabilities.
