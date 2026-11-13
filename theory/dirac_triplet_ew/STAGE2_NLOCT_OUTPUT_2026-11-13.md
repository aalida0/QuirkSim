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

**STATUS: NLOCT_OUTPUT_GENERATED / CONTENT_AUDIT_PENDING /
NLO_EW_PHYSICS_NOT_CERTIFIED.**
