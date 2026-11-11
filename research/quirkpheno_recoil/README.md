# Experimental public-backend weak-recoil pilot

This adapter is a patch against tvami/QuirkPheno commit
`d91cb27ed3786bc63900faad892513841b5b6ed4`. Upstream is Apache-2.0;
its NOTICE credits the ATLAS quirk extension and Geant4 transportation.
Keep upstream LICENSE/NOTICE in the checked-out source. This directory stores
an integration patch and new source, not a second independent transport solver.

**Status: native core and synthetic audit tested; full Geant4 compilation and
transport/causal validation pending AF execution.** Do not merge into stable code.

The patch registers neutral Quirk definitions with sign-preserving endpoint
slots, transfers parent ownership through a one-use token, suppresses ordinary
endpoint death during that transfer, and bypasses production reinitialization
on the neutral daughter's first start. Fingerprints check both incoming string
histories, initialization/borrowed state and liveness across the handoff.
Cached endpoint velocity and kinetic energy are updated to the daughter.
No string restart, reset, projection or partner kick is introduced.

A post-step process performs a physical charged-parent -> neutral-daughter +
pion decay at the first accepted step at/above proper time 50 ps. This is a
**fixed-threshold pion integration test**, not a stochastic lifetime ensemble.
Threshold overshoot and the upstream approximate proper-clock integration need
resolution review. The parent neutral mass is 225 GeV; the parent charged
mass is 225.1571702763212 GeV from the existing QuirkSim LO triplet calculation.
The pion mass used in Geant4 is read from its particle table and logged.

The two rest-frame pion directions +/-y give physical angular controls. Neither
is a no-recoil limit. A stable-versus-decaying comparison mixes charge, mass,
and recoil. This pilot does not implement leptonic channels or all multiplet
production modes and does not certify the Dirac-triplet model globally.

## AF runner

`python3 scripts/af_quirkpheno_weak_pilot.py --output-dir <new-directory>`

Uses an isolated fresh upstream checkout under the output directory and applies
this patch plus the added headers/source. The earlier stable checkout is not
modified. Manifest records upstream revision, adapter file hashes, macros and
execution results. LCG_110, Geant4 11.4.1, one CPU, B=2 T, Lambda=774 eV,
point production, 130 ps. Eight cases at levels 7/8:
- stable charged parents;
- endpoint-0 pion +y and -y;
- endpoint-1 pion +y (charge-conjugate endpoint transfer).

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
