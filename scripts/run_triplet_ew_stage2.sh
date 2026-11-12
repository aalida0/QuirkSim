#!/usr/bin/env bash
# Run full QCD+EW NLOCT with fail-closed log and output checks.
# Does not certify an NLO-EW physics model.
set -Eeuo pipefail
cd "$(dirname "$0")/.."

LOG="${1:-results/triplet_ew_stage2.log}"
NLO="theory/dirac_triplet_ew/generated/DiracTripletEW_QCDEW_CT.nlo"
AUDIT="theory/dirac_triplet_ew/generated/stage2_log_audit.json"
mkdir -p "$(dirname "$LOG")"
rm -f "$AUDIT"
if [[ -e "$NLO" ]]; then
  echo "ERROR: existing NLOCT output $NLO; preserve it before another run." >&2
  exit 2
fi
if [[ -e "$LOG" ]]; then
  previous_log="${LOG}.previous_$(date -u +%Y%m%dT%H%M%SZ)_${BASHPID}"
  mv -- "$LOG" "$previous_log"
  echo "Previous NLOCT stage-2 log preserved: $previous_log"
fi

set +e
wolframscript -file theory/dirac_triplet_ew/02_calculate_nloct.wls \
  2>&1 | tee "$LOG"
run_status=${PIPESTATUS[0]}
set -e

set +e
python3 scripts/audit_triplet_ew_stage2.py \
  --log "$LOG" --nlo "$NLO" --output "$AUDIT"
audit_status=$?
set -e
if ((run_status != 0 || audit_status != 0)); then
  rm -f "$AUDIT"
  echo "NLOCT_STAGE2_REJECTED: Mathematica=$run_status audit=$audit_status" >&2
  echo "Do NOT run stage 03 until the errors and underlying physics are understood." >&2
  exit 1
fi
echo "NLOCT_STAGE2_PREFLIGHT_PASS_ONLY: counterterms require independent EW validation."
