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
