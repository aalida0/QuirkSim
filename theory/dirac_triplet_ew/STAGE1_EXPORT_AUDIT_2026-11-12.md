# QuirkSim NLO EW FeynRules stage-1 failure (2026-11-12)

**Outcome: BLOCKED — NOT NLOCT-READY.** User's FeynRules 2.3.49 /
FeynArts 3.11 run generated `DiracTripletEW_QCDEW_FA.mod`
and `DiracTripletEW_QCDEW_FA.pars`, but logged errors before writing.

## Evidence from user stage1 log
- Merged SM+Triplet model loaded.
- Triplet hermiticity test passed.
- `OnShellRenormalization[...,QCDOnly->False]` produced
  `Part::partd` on `Null` during tadpole redefinition.
- FeynArts writer encountered `Function::slotn`,
  `StringJoin::string`, `Part::partw`, `Set::partd`
  and `First::nofirst`; explicit ignored unsupported
  `V` Lorentz vertices.
- Warned `LorentzTadpole.gen` required for one-point
  scalar vertices, and reported single-W charge nonconservation.
- Printed `FR_STAGE1_FEYNARTS_MODEL_OK` based solely
  on .mod existence — a now-corrected false positive.

The output is NOT suitable for an NLOCT EW calculation.
The presence of a .mod file does not demonstrate that its
UV/R2 vertices, tadpole terms or EW charge renormalization
are correct.

## Remediations pushed to feature/dirac-triplet-ew-lo
- `scripts/audit_triplet_ew_stage1.py`: fail-closed
  Mathematica log and saved model SHA256 screen.
- `scripts/run_triplet_ew_stage1.sh`: archives the
  old output before attempting another stage1;
  audits the log and refuses to leave a pass manifest
  if Mathematica failed or emitted broken vertices.
- `theory/dirac_triplet_ew/01_prepare_feynarts.wls`:
  explicitly stops when EW renormalization emits invalid
  Mathematica expressions.
- `theory/dirac_triplet_ew/02_calculate_nloct.wls`:
  requires a passing manifest with matching model SHA256.
- `theory/dirac_triplet_ew/00_tree_vertex_probe.wls`:
  isolates the tree-level triplet gauge interactions
  from counterterm/SM tadpole generation.

## Next: isolate the triplet vertex model
```bash
cd ~/QuirkSim
git pull --ff-only origin feature/dirac-triplet-ew-lo
export FEYNRULES_HOME="$HOME/.Wolfram/Applications/FeynRules"
set -o pipefail
wolframscript -file theory/dirac_triplet_ew/00_tree_vertex_probe.wls \
  2>&1 | tee results/triplet_ew_tree_probe.log
grep -E 'FR_TREE|FR_TRIPLET|ERROR|::|not conserved' \
  results/triplet_ew_tree_probe.log | tail -80
```

Only if the tree-level vertices agree with the validated
QCD-only UFO should the EW tadpole, generic Lorentz,
wavefunction and parameter renormalization be repaired.
No NLO-EW cross section can yet be quoted.

See also official FeynRules FeynArts interface documentation:
https://cp3.irmp.ucl.ac.be/projects/feynrules/wiki/FeynArtsInterface
and FeynRules' `LorentzTadpole.gen`:
https://github.com/HEPcodes/FeynRules/blob/master/Interfaces/FeynArts/LorentzTadpole.gen


## Follow-up: triplet-only FeynmanRules returned zero vertices

After the stage-1 file-writing errors, the user ran
`00_tree_vertex_probe.wls` and obtained:

```text
FR_TRIPLET_CLASS_DECLARATIONS_PRESENT: {True,True,True}
FR_TREE_VERTEX_COUNT: 0
ERROR: fewer than 8 triplet tree vertices; cannot proceed to EW NLO.
```

A **concrete source error** was identified in
`TripletDirac.fr`: the Dirac fermion conjugates had been
written as `bar[Qp1]`, `bar[Qzero]`, `bar[Qm2]`,
which are NOT the normal FeynRules antiparticle field
symbols. FeynRules creates `Qp1bar`, `Qzerobar` and
`Qm2bar` for non-self-conjugate fermions. (The supported
function alternative is `anti[Qp1]`.) The source now uses
the proper generated conjugates throughout the kinetic,
mass and eight gauge-interaction terms.

The reason for zero Feynman rules has a **strong candidate
explanation**, but its resolution remains unconfirmed
until Mathematica re-runs. The tree probe now reports the
evaluated `LTripletDirac` expression preview and verifies
that it is nonzero, before calling `FeynmanRules`. The
stage-1 renormalization driver also refuses to proceed
unless it finds at least eight tree-level triplet vertices.

Regression test:
`tests/test_triplet_feynrules_antifield_syntax.py`.

**Required next run:**
```bash
cd ~/QuirkSim
git pull --ff-only origin feature/dirac-triplet-ew-lo
export FEYNRULES_HOME="$HOME/.Wolfram/Applications/FeynRules"
set -o pipefail
wolframscript -file theory/dirac_triplet_ew/00_tree_vertex_probe.wls \
  2>&1 | tee results/triplet_ew_tree_probe.log
grep -E 'FR_TREE|FR_TRIPLET|ERROR|::|not conserved' \
  results/triplet_ew_tree_probe.log | tail -80
```

Do not proceed to stage 01 full EW-renormalization or
stage 02 NLOCT until the tree interaction vertices have
been independently validated.

Reference: FeynRules 2.0 paper, arXiv:1310.1921, Section
on fermion class naming and anti[psi].

## Follow-up: eight triplet tree vertices PASS (user Mathematica output)

The corrected `TripletDirac.fr` with proper `Qp1bar`,
`Qzerobar`, `Qm2bar` generated exactly **8**
`FeynmanRules[LTripletDirac,MaxParticles->3]` FFV vertices,
consistent with the previously validated NLO QCD model:

| Fermion | Antifermion | Boson | Coefficient / gamma^mu |
|---|---|---|---|
| Qm2 | Qm2bar | A | `-i ee` |
| Qp1 | Qp1bar | A | `+i ee` |
| Qp1 | Qzerobar | Wbar | `+i ee/sw` |
| Qzero | Qp1bar | W | `+i ee/sw` |
| Qm2 | Qzerobar | W | `+i ee/sw` |
| Qzero | Qm2bar | Wbar | `+i ee/sw` |
| Qm2 | Qm2bar | Z | `-i cw ee/sw` |
| Qp1 | Qp1bar | Z | `+i cw ee/sw` |

The tree FFV test has **no evidence of a missing or unexpected
coupling**. It does not independently establish gauge Ward identities
or full one-loop EW renormalization.

Additional implemented checks:
- `scripts/audit_triplet_ew_tree_vertices.py` parses the Mathematica
  log and checks the eight FFV coefficients against hard-coded validated
  QCD-UFO expectations, refusing missing, duplicated, wrong-sign,
  chiral or wrong-normalization vertices;
- `scripts/run_triplet_ew_stage1.sh` now runs that tree audit and
  aborts before renormalization if it fails.
- `01_prepare_feynarts.wls` now uses `GenericFile -> True`, generating
  a **matching model-specific .gen** rather than filtering every
  interaction against `Lorentz.gen`, which previously ignored vertices.
- `02_calculate_nloct.wls` now uses the corresponding model-specific
  generic file, and the stage-1 preflight requires .mod, .pars and
  .gen with matching SHA256 hashes.

**Next local action**: after pulling the branch, run
`bash scripts/run_triplet_ew_stage1.sh` (NOT stage2). It
re-runs the verified tree probe and attempts EW renormalization
with a fail-closed log audit. The prior defective EW stage1 export
was never physics-validated, and this new run may expose further
SM tadpole renormalization issues.


## Follow-up: stage-1 model-specific FeynArts export PASSES the guarded preflight

Latest user-provided stage1 driver output:
- `FeynRules 2.3.49` with `FeynArts 3.11`.
- Tree FFV audit: `TREE_FFV_UFO_COMPATIBLE_PRELIMINARY`,
  eight vertices and no unexpected couplings.
- `FR_STAGE1_TRIPLET_TREE_VERTEX_COUNT: 8`.
- `CheckHermiticity`: **hermitian**.
- `OnShellRenormalization[...,QCDOnly->False]`: this run
  had no previously observed fatal Mathematica messages. It
  emitted `Solve::nongen` warnings (nonfatal generic solution warnings),
  and completed the tadpole redefinition.
- FeynArts writer: **203 generated vertices** from **219**
  possible candidate structures, with no previously seen
  unsupported Lorentz-gen or ignored-vertex messages.
- Files: `DiracTripletEW_QCDEW_FA.{mod,pars,gen}` written
  under `theory/dirac_triplet_ew/generated/`.
- Stage1 audit: `STAGE1_LOG_PREFLIGHT_PASS_NOT_FULL_PHYSICS_CERTIFIED`,
  all output files present, `blocking_findings: {}`.
- Status: `STAGE1_LOG_PREFLIGHT_OK`. This is a **source and
  structural preflight**, not validation of loop EW UV/R2.

Validated generated model SHA256s (for the user-provided run):
- .mod `4783c579beb2ffcbbda282f50481fbf43f6aa27fd1528c787b16e1941fb28e6f`
- .pars `04914d33dbb21b471a0b816abc133cc3495d9a2f56110835adb425763940e3a3`
- .gen `ecebe86804650719d9dc733a1232eec875cf8e746639078fd734e7311721dc6f`

### Move to NLOCT stage 02

The prepared driver `scripts/run_triplet_ew_stage2.sh` now:
- refuses to overwrite earlier stage-2 logs or `.nlo` files,
- runs NLOCT in a fresh Wolfram kernel with
  `WriteCT[classModel, matchingGeneric, QCDOnly->False]`,
- verifies a successful stage1 log audit and exact SHA256 of
  `.mod`, `.pars`, `.gen`,
- checks Mathematica structural diagnostics, and verifies
  the UV and R2 lists and model provenance in the `.nlo`,
- emits a stage2 JSON gate used by stage 03.

Commands:

```bash
cd ~/QuirkSim
git pull --ff-only origin feature/dirac-triplet-ew-lo
export FEYNRULES_HOME="$HOME/.Wolfram/Applications/FeynRules"
export FEYNARTS_HOME="$HOME/.Wolfram/Applications/FeynArts-3.11"
bash scripts/run_triplet_ew_stage2.sh
grep -E 'NLOCT_STAGE2|ERROR|::|Writing R2 and UV' \
  results/triplet_ew_stage2.log | tail -80
```

**Important:** The NLOCT stage is potentially computationally demanding
because `QCDOnly->False` requests EW plus QCD loop sectors of the
complete SM and BSM model. The NLOCT engine and renormalization
scheme need physical scrutiny after any numerical/export success.
Even a stage2 preflight pass does NOT validate
`NLOEW` cross sections or authorize publication.


## NLOCT Stage 2 initial run: FeynArts model lookup failure

The first stage-2 attempt ran FeynArts 3.11 and NLOCT 1.02,
passed the saved stage-1 SHA256 check, then failed immediately at
`Computing CT for the {S} vertices` with

```text
$Aborted
InitializeModel::noopen: Cannot open DiracTripletEW_QCDEW_FA/DiracTripletEW_QCDEW_FA.gen.
NLOCT_STAGE2_REJECTED
nlo_size_bytes: null
```

**Root cause:** the generated class/generic model is stored in a nested
directory under `theory/dirac_triplet_ew/generated/`. FeynArts
`ReadModelFile` searches `Block[{$Path = $ModelPath},
System`Private`FindFile[model <> ".gen"]]`, so
`SetDirectory[generated]` alone does not ensure model lookup.
The .gen hash is verified and this is a **path resolution failure**,
not evidence about UV/R2 physics.

**Fix:** In `02_calculate_nloct.wls`, prepend the absolute generated
model subdirectory to FeynArts `$ModelPath`, resolve both
`DiracTripletEW_QCDEW_FA.gen` and `.mod` using FeynArts'
own `FindFile` convention, explicitly check
`InitializeModel[qsModelName, GenericModel -> qsModelName]`,
and pass just that common basename to `WriteCT`. The script
aborts before the heavy CT calculation if lookup or initialization
fails. The shell driver `run_triplet_ew_stage2.sh` automatically
archives the failed previous stage2 log before retrying, while
refusing to overwrite any actual `.nlo` output.

**Next:** pull `feature/dirac-triplet-ew-lo` and rerun
`bash scripts/run_triplet_ew_stage2.sh` with the existing
`FEYNRULES_HOME` and `FEYNARTS_HOME` exports. Confirm
`NLOCT_GENERIC_MODEL_RESOLVED` and
`NLOCT_FEYNARTS_INITIALIZATION_OK` before the slow UV/R2
calculation starts. Only NLOCT .nlo generation and the subsequent
UV/IR gauge tests can establish physical readiness.

Primary source: FeynArts `Initialize.m` source, `ReadModelFile`
and `$ModelPath`:
https://github.com/FeynCalc/feynarts-mirror/blob/master/FeynArts/Initialize.m
