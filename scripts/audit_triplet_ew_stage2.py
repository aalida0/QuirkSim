#!/usr/bin/env python3
"""Fail-closed preflight for FeynArts/NLOCT stage-2 QCD+EW .nlo output.

This checks only the generation log and basic presence/syntax of two
counterterm collections. It NEVER certifies the UV/R2 physics, gauge
or QED-IR cancellation, hidden-color multiplicities or EW scheme.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path

FATAL = {
    "mathematica_structure_error":
      r"\b(?:Part|Set)::part(?:d|w)\b|\b(?:First::nofirst|Last::nolast|Function::slotn|StringJoin::string)\b",
    "nonfinite_expression":r"\b(?:Power::infy|Infinity::indet)\b|\bIndeterminate\b",
    "mathematica_abort":r"\$Aborted|\bAbort\[",
    "explicit_error":r"\bERROR:|\bNLOCT_STAGE2_REJECTED",
    "ignored_vertex":r"Vertex is ignored\.|not supported by Lorentz(?:Tadpole)?\.gen",
}
SIGNALS=(
    "FR_STAGE1_LOG_AUDIT_PASSED_ONLY",
    "NLOCT_STAGE2_FILE_WRITTEN_UNVALIDATED",
)


def audit(log_path: Path,nlo_path: Path) -> dict:
    log_path=log_path.resolve()
    nlo_path=nlo_path.resolve()
    raw=log_path.read_bytes()
    log=raw.decode("utf-8",errors="replace")
    findings=[]
    for label,regex in FATAL.items():
        if re.search(regex,log,re.I|re.S):
            findings.append(f"Rejected Mathematica/NLOCT log pattern: {label}")
    for sig in SIGNALS:
        if sig not in log:
            findings.append(f"Missing stage2 success marker: {sig}")
    content=""
    if not nlo_path.is_file():
        findings.append("No .nlo file written")
    else:
        if nlo_path.stat().st_size<200:
            findings.append("NLOCT .nlo file is too small for a renormalized EW model")
        content=nlo_path.read_text(encoding="utf-8",errors="replace")
        for var in ("UV$vertlist","R2$vertlist"):
            match=re.search(re.escape(var)+r"\s*=\s*\{",content)
            if not match:
                findings.append(f"Missing counterterm list: {var}")
            elif re.search(re.escape(var)+r"\s*=\s*\{\s*\}\s*;",content):
                findings.append(f"Empty counterterm list: {var}")
        if "FR$InteractionOrderPerturbativeExpansion" not in content:
            findings.append("No NLOCT perturbative-order metadata")
        if "CT$Model" not in content or "CT$GenericModel" not in content:
            findings.append("Missing provenance for FeynArts class/generic model")
        if re.search(r"\b(?:Indeterminate|ComplexInfinity)\b",content):
            findings.append("Nonfinite expression found in NLOCT .nlo")
    return {
        "status":"NLOCT_STAGE2_STRUCTURAL_PREFLIGHT_PASS_NOT_PHYSICS_CERTIFIED"
            if not findings else "NLOCT_STAGE2_REJECTED",
        "preflight_pass":not findings,
        "log":str(log_path),
        "log_sha256":hashlib.sha256(raw).hexdigest(),
        "nlo":str(nlo_path),
        "nlo_size_bytes":nlo_path.stat().st_size if nlo_path.is_file() else None,
        "nlo_sha256":hashlib.sha256(nlo_path.read_bytes()).hexdigest()
            if nlo_path.is_file() else None,
        "blocking_findings":findings,
        "warning":"Passing only confirms plausible output; EW scheme, UV/R2 coverage and cancellation unverified.",
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log",type=Path,required=True)
    parser.add_argument("--nlo",type=Path,
        default=Path("theory/dirac_triplet_ew/generated/DiracTripletEW_QCDEW_CT.nlo"))
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=audit(args.log,args.nlo)
    print(json.dumps(result,indent=2))
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    if not result["preflight_pass"]:
        raise SystemExit("STOP: NLOCT stage2 output is incomplete or reports Mathematica errors.")
    print("NLOCT_STAGE2_PREFLIGHT_PASS_ONLY: independent UV/IR physics validation remains required.")


if __name__=="__main__":
    main()
