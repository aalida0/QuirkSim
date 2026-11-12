#!/usr/bin/env python3
"""Fail-closed audit of tree-level FeynRules Dirac-triplet FFV vertices.

Checks a saved output from 00_tree_vertex_probe.wls against the
eight FFV vertices and couplings independently implemented in the
validated DiracTripletEW_NLO_UFO source.

This is a tree-level algebraic regression check, not a demonstration
of full SU(2) Ward identities or NLO EW renormalization.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path

# (antifermion, fermion, gauge field): sign, normalization
# For the minimal Y=0 Dirac (Q1+, Q0, Q2-) triplet.
EXPECTED = {
    ("Qm2bar", "Qm2", "A"): (-1, "ee"),
    ("Qp1bar", "Qp1", "A"): (+1, "ee"),
    ("Qp1bar", "Qzero", "W"): (+1, "ee/sw"),
    ("Qzerobar", "Qm2", "W"): (+1, "ee/sw"),
    ("Qm2bar", "Qzero", "Wbar"): (+1, "ee/sw"),
    ("Qzerobar", "Qp1", "Wbar"): (+1, "ee/sw"),
    ("Qm2bar", "Qm2", "Z"): (-1, "cw*ee/sw"),
    ("Qp1bar", "Qp1", "Z"): (+1, "cw*ee/sw"),
}
VERTEX_RE = re.compile(
    r"FR_TREE_VERTEX_(\d+):\s*InputForm\["
    r"\{\{\{(\w+),\s*1\},\s*\{(\w+),\s*2\},\s*\{(\w+),\s*3\}\},\s*"
    r"(.+)\}\]\s*$",
    flags=re.MULTILINE,
)
FATAL_RE = re.compile(r"\bERROR:|\$Aborted|\b(?:Part|Set|First|Function|StringJoin)::",re.I)


def _coeff_ok(expr: str, sign: int, kind: str) -> bool:
    # The Wolfram InputForm representation of the FFV vertex is unambiguous:
    # coefficient (I*ee, I*ee/sw, I*cw*ee/sw) times one gamma tensor.
    src="".join(expr.split())
    if src.count("Ga[")!=1:
        return False
    if "Index[Lorentz,Ext[3]]" not in src:
        return False
    if "Index[Spin,Ext[1]]" not in src or "Index[Spin,Ext[2]]" not in src:
        return False
    if re.search(r"Proj[MP]|Gamma5|PL|PR",src):
        return False
    tensor="Ga[Index[Lorentz,Ext[3]],Index[Spin,Ext[1]],Index[Spin,Ext[2]]]"
    # All variations previously observed (parentheses optional)
    if kind=="ee":
        allowed=("-I*ee*"+tensor, "I*ee*"+tensor)
    elif kind=="ee/sw":
        allowed=("(-I*ee*"+tensor+")/sw","(I*ee*"+tensor+")/sw")
    elif kind=="cw*ee/sw":
        allowed=("(-I*cw*ee*"+tensor+")/sw","(I*cw*ee*"+tensor+")/sw")
    else:
        raise ValueError("unexpected expected factor")
    return src==allowed[1 if sign>0 else 0]


def audit(log_path: Path) -> dict:
    log_path=log_path.resolve()
    raw=log_path.read_bytes()
    data=raw.decode("utf-8",errors="replace")
    seen={}
    findings=[]
    head_signals=(
        "FR_TRIPLET_CLASS_DECLARATIONS_PRESENT:",
        "FR_TRIPLET_LAGRANGIAN_DECLARATION_OK",
        "FR_TREE_PROBE_START",
        "FR_TREE_PROBE_COMPLETED_UNVALIDATED",
    )
    for marker in head_signals:
        if marker not in data:
            findings.append("Missing expected success marker: "+marker)
    declared=re.search(r"FR_TREE_VERTEX_COUNT:\s*(\d+)",data)
    if declared is None or int(declared.group(1))!=8:
        findings.append("FeynRules did not report exactly eight triplet tree vertices")
    if FATAL_RE.search(data):
        findings.append("Wolfram error/abort detected in probe log")
    found=list(VERTEX_RE.finditer(data))
    if len(found)!=8:
        findings.append(f"Expected exactly 8 parsed vertex records, got {len(found)}")
    for m in found:
        index=int(m.group(1))
        key=tuple(m.group(i) for i in (2,3,4))
        coeff=m.group(5)
        if key in seen:
            findings.append(f"Duplicate process vertex: {key}")
        seen[key]=index
        if key not in EXPECTED:
            findings.append(f"Unexpected vertex: {key}")
            continue
        sign,factor=EXPECTED[key]
        if not _coeff_ok(coeff,sign,factor):
            findings.append(f"Unexpected coefficient/tensor in vertex {index}: {key}: {coeff}")
    missing=set(EXPECTED)-set(seen)
    if missing:
        findings.append("Missing vertices: "+str(sorted(missing)))
    return {
        "status": "TREE_FFV_UFO_COMPATIBLE_PRELIMINARY" if not findings
            else "TREE_FFV_AUDIT_REJECTED",
        "passed": not findings,
        "log":str(log_path),
        "log_sha256":hashlib.sha256(raw).hexdigest(),
        "declared_vertices":int(declared.group(1)) if declared else None,
        "parsed_vertices":len(found),
        "missing_or_unexpected":findings,
        "warning":"Check couplings against source conventions. This does not validate EW loop counterterms.",
    }


def main() -> None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--log",type=Path,required=True)
    ap.add_argument("--output",type=Path)
    a=ap.parse_args()
    result=audit(a.log)
    print(json.dumps(result,indent=2))
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    if not result["passed"]:
        raise SystemExit("REJECTED: triplet tree vertices differ from validated QCD UFO")
    print("PASS: eight Dirac-triplet FFV vertices match the NLO-QCD UFO (tree level).")


if __name__=="__main__":
    main()
