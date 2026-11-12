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
