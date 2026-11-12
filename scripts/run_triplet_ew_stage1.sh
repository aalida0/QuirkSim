#!/usr/bin/env bash
# Run FeynRules stage1 and fail closed if Mathematica emitted hidden errors.
set -Eeuo pipefail

cd "$(dirname "$0")/.."
mkdir -p results
TREE_LOG="results/triplet_ew_tree_probe.log"
TREE_AUDIT="theory/dirac_triplet_ew/generated/tree_vertex_audit.json"
LOG="${1:-results/triplet_ew_stage1.log}"
MODEL="theory/dirac_triplet_ew/generated/DiracTripletEW_QCDEW_FA"
AUDIT="theory/dirac_triplet_ew/generated/stage1_log_audit.json"
mkdir -p "$(dirname "$LOG")"
rm -f "$TREE_AUDIT" "$AUDIT"

# Hard tree-level validation first; failures must never invoke EW
# OnShellRenormalization or overwrite a previous FeynArts model.
set +e
wolframscript -file theory/dirac_triplet_ew/00_tree_vertex_probe.wls \
  2>&1 | tee "$TREE_LOG"
tree_status=${PIPESTATUS[0]}
python3 scripts/audit_triplet_ew_tree_vertices.py \
  --log "$TREE_LOG" --output "$TREE_AUDIT"
tree_audit_status=$?
set -e
if ((tree_status != 0 || tree_audit_status != 0)); then
  rm -f "$TREE_AUDIT"
  echo "TREE_FFV_REJECTED: Wolfram=$tree_status; audit=$tree_audit_status" >&2
  echo "Do NOT attempt the electroweak renormalization." >&2
  exit 1
fi
# Preserve any earlier FeynArts output so a stale .mod cannot pass this run.
if [[ -d "$MODEL" ]]; then
  old_name="${MODEL}.previous_$(date -u +%Y%m%dT%H%M%SZ)_$"
  mv "$MODEL" "$old_name"
  echo "Previous FeynArts model preserved: $old_name"
fi
# Never leave a stale previous PASS if this run is interrupted or fails.
rm -f "$AUDIT"

set +e
wolframscript -file theory/dirac_triplet_ew/01_prepare_feynarts.wls \
  2>&1 | tee "$LOG"
run_status=${PIPESTATUS[0]}
set -e

# Always write a log audit JSON (even if Wolfram exits nonzero). This fails on
# errors/ignored vertices that do not affect Wolfram's process exit code.
set +e
python3 scripts/audit_triplet_ew_stage1.py \
  --log "$LOG" --model-dir "$MODEL" --json "$AUDIT"
audit_status=$?
set -e
if ((run_status != 0 || audit_status != 0)); then
  rm -f "$AUDIT"  # reject even if Mathematica failed but static log was accidentally clean
  echo "STAGE1_REJECTED: Mathematica exit status=$run_status; log audit exit status=$audit_status" >&2
  echo "Do NOT proceed to NLOCT stage 02." >&2
  exit 1
fi

echo "STAGE1_LOG_PREFLIGHT_OK only. Verify the physics before treating the model as UV/R2-complete."
