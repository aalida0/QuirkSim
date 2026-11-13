#!/usr/bin/env bash
# Guarded candidate UFO export. Do NOT use candidate for physics without
# NLO-EW UV/IR cancellation, Ward tests, SM limit, and N_IC loop audit.
set -Eeuo pipefail
cd "$(dirname "$0")/.."
mkdir -p results
WOLFRAMSCRIPT_BIN="${WOLFRAMSCRIPT_BIN:-wolframscript}"
if ! command -v "$WOLFRAMSCRIPT_BIN" >/dev/null 2>&1; then
  echo "ERROR: no wolframscript found at WOLFRAMSCRIPT_BIN=$WOLFRAMSCRIPT_BIN" >&2
  exit 2
fi

LOG="${1:-results/triplet_ew_stage3.log}"
UFO="theory/dirac_triplet_ew/generated/DiracTripletEW_QCDEW_CANDIDATE_UFO"
if [[ -e "$UFO" ]]; then
  echo "ERROR: candidate UFO already exists: $UFO" >&2
  echo "Archive it explicitly before another export; will not overwrite." >&2
  exit 2
fi
if [[ -f "$LOG" ]]; then
  earlier="${LOG}.previous_$(date -u +%Y%m%dT%H%M%SZ)_${BASHPID}"
  mv -- "$LOG" "$earlier"
  echo "Prior stage-3 log saved: $earlier"
fi

echo "WARNING: importing a 136-MiB NLOCT file may require substantially"
echo "more RAM in Mathematica than its on-disk size. Never treat output"
echo "presence as NLO-EW physics validation."
set +e
"$WOLFRAMSCRIPT_BIN" -file theory/dirac_triplet_ew/03_export_candidate_ufo.wls \
  2>&1 | tee "$LOG"
rc=${PIPESTATUS[0]}
set -e
if ((rc != 0)); then
  echo "NLOCT_STAGE3_REJECTED: Mathematica exit status $rc" >&2
  exit 1
fi
if grep -Eq '(NLOCT_STAGE3_REJECTED|ERROR:|\\$Aborted|InitializeModel::noopen|SetDirectory::badfile|SetDirectory::cdir|Vertex is ignored)' "$LOG"; then
  echo "NLOCT_STAGE3_REJECTED: see log for missing files/Mathematica errors" >&2
  exit 1
fi
if ! grep -q 'NLOCT_STAGE3_CANDIDATE_FILES_PRESENT_UNVALIDATED' "$LOG"; then
  echo "NLOCT_STAGE3_REJECTED: success marker absent" >&2
  exit 1
fi

echo "CANDIDATE_UFO_EXPORTED_UNVALIDATED at $UFO"
echo "DO NOT quote NLO EW rates before MG5 UV/IR/gauge tests."
