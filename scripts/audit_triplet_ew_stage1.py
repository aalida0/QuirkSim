#!/usr/bin/env python3
"""Fail-closed audit of the FeynRules -> FeynArts full-EW stage-1 log.

FeynRules may write a .mod/.pars after silently discarding vertices, or after
non-fatal Mathematica messages in its loop renormalization. Never treat file
existence alone as permission to invoke NLOCT (stage 2).

This log screen does not certify a valid physics model: a PASS only means no
recognized fatal messages were encountered. Independent vertices and the
EW-renormalization conditions still need review before stage 2 results
can be used for cross sections.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

# Unresolved structural/Mathematica errors that invalidate this FeynArts export.
FATAL_PATTERNS = {
    "mathematica_part_error": r"\b(?:Part|Set)::part(?:d|w)\b",
    "mathematica_slot_error": r"\bFunction::slotn\b",
    "mathematica_string_error": r"\bStringJoin::string\b",
    "mathematica_first_error": r"\bFirst::nofirst\b",
    "counterterm_renormalization_tadpole_error": r"Tadpole redefinition[\s\S]{0,300}?Part::partd",
    "unrecognized_generic_coupling": r"coupling .*?not supported by Lorentz(?:Tadpole)?\.gen",
    "ignored_vertex": r"Vertex is ignored\.",
    "unsupported_tadpole": r"tadpole vertices .*?this model should be used",
    "charge_nonconservation": r"Quantum number Q not conserved in vertex",
    "aborted": r"(?:\$Aborted|\bAbort\[|\bERROR:)",
}
PASS_SIGNALS = (
    "FR_TRIPLET_MERGED_MODEL_LOADED",
    "The lagrangian is hermitian.",
)


def audit(log_path: Path, model_dir: Path) -> dict:
    log_path = log_path.resolve()
    model_dir = model_dir.resolve()
    text = log_path.read_text(encoding="utf-8",errors="replace")
    errors = {}
    for name,pattern in FATAL_PATTERNS.items():
        # regexp multi-line when a Mathemtica printer wrapped text across lines
        matched = list(re.finditer(pattern,text,re.IGNORECASE|re.DOTALL))
        if matched:
            errors[name]={
                "count": len(matched),
                "samples": [text[max(0,m.start()-50):m.end()+50].replace("\n"," ")[:180]
                            for m in matched[:3]],
            }
    # Recent stage 1 uses a "FILES_WRITTEN" banner; accept historical 'OK' banner
    # only as proof that the files were written, not as proof of correctness.
    exists={s:(model_dir/(model_dir.name+s)).is_file()
            for s in (".mod",".pars")}
    signals={s:(s in text) for s in PASS_SIGNALS}
    signals["export_reported"]=(
        "FR_STAGE1_FEYNARTS_MODEL_OK" in text or
        "FR_STAGE1_FEYNARTS_FILES_WRITTEN_UNVALIDATED" in text)
    accepted=(not errors and all(exists.values()) and all(signals.values()))
    return {
        "status": "STAGE1_LOG_PREFLIGHT_PASS_NOT_FULL_PHYSICS_CERTIFIED" if accepted
                  else "STAGE1_REJECTED_UNSAFE_FOR_NLOCT",
        "accepted_for_nloct_preflight": accepted,
        "log": str(log_path),
        "model_dir":str(model_dir),
        "log_sha256": hashlib.sha256(log_path.read_bytes()).hexdigest(),
        "model_sha256": {s:hashlib.sha256(
            (model_dir/(model_dir.name+s)).read_bytes()).hexdigest()
            for s,present in exists.items() if present},
        "signals":signals,
        "model_files_exist":exists,
        "blocking_findings":errors,
        "warning":"A passing log audit does not establish complete one-loop UV/R2 physics.",
    }


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--log",type=Path,required=True)
    p.add_argument("--model-dir",type=Path,
        default=Path("theory/dirac_triplet_ew/generated/"
                     "DiracTripletEW_QCDEW_FA"))
    p.add_argument("--json",type=Path)
    a=p.parse_args()
    report=audit(a.log,a.model_dir)
    print(json.dumps(report,indent=2))
    if a.json:
        a.json.parent.mkdir(parents=True, exist_ok=True)
        a.json.write_text(json.dumps(report,indent=2)+"\n")
    if not report["accepted_for_nloct_preflight"]:
        raise SystemExit("STOP: FeynArts files were written, but stage-1 contains"
                         " invalid/ignored vertices or lacks mandatory markers."
                         " Do not run NLOCT.")
    print("PASS: static stage-1 log preflight only; loop physics still unvalidated.")


if __name__=="__main__":
    main()
