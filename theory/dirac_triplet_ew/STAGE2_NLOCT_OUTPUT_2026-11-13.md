# NLOCT full QCD+EW stage 2 output: generated but not physics validated

**Date:** 2026-11-13. **Branch:** `feature/dirac-triplet-ew-lo`.

## Evidence from the workstation output

The user ran NLOCT 1.02 with FeynArts 3.11 on the generated
FeynRules SM+Dirac-triplet full-counterterm FeynArts model.
The previous model-lookup problem has been resolved:

```text
NLOCT_GENERIC_MODEL_RESOLVED: .../DiracTripletEW_QCDEW_FA.gen
NLOCT_CLASSES_MODEL_RESOLVED: .../DiracTripletEW_QCDEW_FA.mod
NLOCT_FEYNARTS_INITIALIZATION_OK
NLOCT_STAGE2_FILE_WRITTEN_UNVALIDATED:
  .../generated/DiracTripletEW_QCDEW_CT.nlo (142907954 bytes)
NLOCT_STAGE2_REQUIRES_LOG_AND_CONTENT_AUDIT
```

**Important**: these markers establish a successfully initialized
FeynArts model and a physically written NLOCT `.nlo` file,
not a demonstration of UV finiteness, gauge invariance,
correct EW input/IR scheme, or NLO-EW cross-section validity.

The user printed selected `results/triplet_ew_stage2.log`
lines, which do NOT include the Python wrapper's separate
`NLOCT_STAGE2_PREFLIGHT_PASS_ONLY` audit result.
Check the saved stage2 audit JSON independently before stage 3.

The file is **142,907,954 bytes (~136 MiB)**. It is far too large
to paste into a chat or commit to regular GitHub (normal GitHub
repository uploads have a 100-MiB hard file-size limit).
The exported UFO and its counterterms may be large, and
the Mathematica step may consume substantially more memory
than the on-disk file; avoid running stage 3 blindly.

## Compact, streaming inspection

A streaming `scripts/inspect_triplet_nloct_large.py` audit
checks file completion, required `R2$vertlist`, `UV$vertlist`
and `FR$CTparam` lists, counts records, scans for nonfinite
tokens and tracks appearances of triplet field names/class IDs
in UV and R2 sectors. It avoids loading the entire 136-MiB
file merely to count entries.

```bash
cd ~/QuirkSim
git pull --ff-only origin feature/dirac-triplet-ew-lo

# First inspect the existing audit JSON from the earlier run:
python3 -m json.tool \
  theory/dirac_triplet_ew/generated/stage2_log_audit.json

# Then run the new streaming full-content inspection:
python3 scripts/inspect_triplet_nloct_large.py \
  --nlo theory/dirac_triplet_ew/generated/DiracTripletEW_QCDEW_CT.nlo \
  --output results/triplet_ew_nloct_content_report.json
```

**Do not run stage 3** until the stage2 audit and the new
content inspector pass, and the reported triplet-specific
R2/UV terms are physically inspected. If the inspector
reports no recognizable triplet fields, that is a manual
investigation flag: FeynArts may remap field indices.
It does not by itself mean the calculated counterterms
are absent.

The stage-3 export script now requires a passing content
report and an exact matching `SHA256` for the `.nlo`
source, in addition to the existing stage2 log audit.

## Physics acceptance checks still required

1. Confirm the UV/R2 lists actually cover all SM + three
   vectorlike triplet fermions in the intended NLO coupling
   orders; no missing tadpole, mixed gauge-kinetic field
   or mass counterterms.
2. Verify UV-pole cancellation for independent
   `u u~ -> q1+ q1-`, `d d~ -> q2- q2+`,
   and charged-neutral channels in the same electroweak
   renormalization scheme and QED factorization convention.
3. Verify genuine QED IR cancellation with real emission,
   photon PDFs, and sensible soft/collinear handling for
   dressed charged final states.
4. Verify weak-mass splitting uses **pole masses** without
   an inadmissible gauge-breaking tree splitting.
5. For `N_IC=2` closed virtual quirk loops carry
   hidden-color trace factors. The present single-copy
   NLOCT calculation is not automatically sufficient.
6. Compare SM-only amplitudes to a validated
   NLO QCD+EW Standard Model; check gauge/Ward identities.
7. Run controlled MG5 fixed-order validation and
   document UV/IR pole cancellations and reproducible
   finite cross sections before quoting NLO QCD+EW.

**STATUS: NLOCT_OUTPUT_GENERATED / TWO_STRUCTURAL_AUDITS_PASS /
EXTERNAL_VERTEX_HEADER_AUDIT_PENDING / NLO_EW_PHYSICS_NOT_CERTIFIED.**


## Verified structural audit from workstation (2026-11-13)

The user ran BOTH audits successfully on the same generated NLOCT file:

```text
NLOCT_STAGE2_STRUCTURAL_PREFLIGHT_PASS_NOT_PHYSICS_CERTIFIED
CONTENT_STRUCTURAL_SCREEN_PASS_NOT_PHYSICS_CERTIFIED
size:      142907954 bytes (~136.29 MiB)
sha256:    58094c1ae9e3127bfdc1bfee4126fe017bd5c4a4cb99e8ec39d22862a9e699ba
lines:     457
R2:        206 approximate record lines
UV:        213 approximate record lines
CTparam:   present; list structurally closed
longest R2 line: 19098 bytes
longest UV line: 2552509 bytes
nonfinite tokens: none
blocking findings: none
```

The `R2` section mentions `Qp1` 61 times, `Qzero` 71,
`Qm2` 61. In the `UV` section, `Qp1` appears 5149 times,
`Qzero` 5127, `Qm2` 4981. **These are whole-expression
token occurrences, NOT independent vertex counts and NOT proof
of triplet external-leg counterterms.** In particular, UV
expressions can mention triplet masses and fields in internal
loops for otherwise SM external legs.

A narrower follow-up tool,
`scripts/inspect_triplet_nloct_vertex_headers.py`, now
streams the `R2$vertlist` and `UV$vertlist` records and
extracts the *external fields* (first item of the NLOCT
`{external-fields, coefficient}` pair) separately from the
large loop coefficient. It reports distinct headers and triplet
external-leg occurrences, with representative headers (never
megabyte expressions). A passing structural parse remains
**NOT a proof of EW gauge or UV/IR consistency**.

Run:

```bash
cd ~/QuirkSim
git pull --ff-only origin feature/dirac-triplet-ew-lo
python3 scripts/inspect_triplet_nloct_vertex_headers.py \
  --nlo theory/dirac_triplet_ew/generated/DiracTripletEW_QCDEW_CT.nlo \
  --output results/triplet_ew_nloct_header_audit.json
```

Stage 03 now requires this matching header-audit manifest,
in addition to both earlier SHA256-matched audit manifests.
**Do not call the candidate UFO an NLO EW prediction** until
amplitude-level UV cancellation, QED real-emission/IR
subtraction, gauge independence, and hidden-index multiplicity
have been demonstrated.

## Candidate UFO export is guarded, not approved for physics

`03_export_candidate_ufo.wls` now requires the stage-2 log audit,
the large-file streaming audit, and the targeted UV/R2 external-leg
header audit to refer to **the exact same .nlo SHA256**. It refuses
overwriting a pre-existing candidate UFO. It checks Mathematica
counterterm loading and the presence of the expected UFO and
`CT_*.py` files, and prints `UNVALIDATED` markers.

An optional guarded shell driver is
`bash scripts/run_triplet_ew_stage3.sh`. The Mathematica
`Get` step needs to materialize the full 136-MiB symbolic NLOCT
input and can consume far more RAM than the file size. No stage-3
run has been reported yet. Do not run it until the header audit
and model-content review have been completed.

The old trusted NLO-QCD UFO is unchanged. No measured
NLO-EW cross section exists.
