#!/usr/bin/env python3
"""Fail-closed static screening of a proposed triplet NLO QCD+EW UFO.

Passing this test only shows the UFO contains some triplet-specific UV/R2
vertices, it does NOT prove that its counterterms or finite mass
renormalization are correct. Publication certification requires the
process-level NLO EW UV/IR tests, gauge/ward identities and independent
cross-checks described in theory/dirac_triplet_ew/README.md.

The existing NLO-QCD-only UFO MUST FAIL this check.
"""
from __future__ import annotations

import argparse
import ast
import json
from collections import Counter
from pathlib import Path


REQUIRED_FILES=(
    "particles.py", "parameters.py", "vertices.py", "couplings.py",
    "coupling_orders.py", "CT_vertices.py", "CT_couplings.py",
    "CT_parameters.py",
)
# FeynRules exports particle class names variously; all three must be found
# in the *tree-level particle file*, and CTs must involve triplet fields.
CHARGED={"q1p","q1m","q2m","q2p","qp1","qm2",
         "q1__plus__","q2__minus__"}
NEUTRAL={"qzero","q0","q0bar","q0__tilde__"}


def _attr_names(node: ast.AST) -> set[str]:
    return {p.attr.lower() for p in ast.walk(node) if isinstance(p, ast.Attribute)}


def inspect_ufo(path: Path) -> dict:
    path=path.resolve()
    missing=[n for n in REQUIRED_FILES if not (path/n).is_file()]
    if missing:
        raise ValueError("Missing mandatory UFO files: "+", ".join(missing))
    particle_tree=ast.parse((path/"particles.py").read_text(encoding="utf-8"))
    ct_tree=ast.parse((path/"CT_vertices.py").read_text(encoding="utf-8"))
    names={n.id.lower() for n in ast.walk(particle_tree) if isinstance(n,ast.Name)}
    if not any(k in names for k in CHARGED) or not any(k in names for k in NEUTRAL):
        raise ValueError("Triplet charged/neutral particle fields missing from UFO")

    counts=Counter()
    details=[]
    for node in ast.walk(ct_tree):
        if not isinstance(node,ast.Call) or not (
            isinstance(node.func,ast.Name) and node.func.id=="CTVertex"
        ):
            continue
        kw={x.arg:x.value for x in node.keywords}
        if "type" not in kw:
            continue
        ct_type=kw["type"]
        if not isinstance(ct_type,ast.Constant) or not isinstance(ct_type.value,str):
            continue
        name=ct_type.value.upper()
        if name not in ("UV","R2","UVMASS","UVLOOP","UVTREE"):
            continue
        parts=_attr_names(kw.get("particles",ast.Constant(value=None)))
        has_charged=bool(parts & CHARGED)
        has_neutral=bool(parts & NEUTRAL)
        if has_charged or has_neutral:
            counts[name]+=1
            details.append({"kind":name,"charged":has_charged,"neutral":has_neutral})
            if has_charged: counts["charged_"+name]+=1
            if has_neutral: counts["neutral_"+name]+=1
    uv_total=sum(v for k,v in counts.items() if k in ("UV","UVMASS","UVLOOP","UVTREE"))
    r2_total=counts["R2"]
    return {"model_path":str(path),"n_triplet_uv_vertices":uv_total,
            "n_triplet_r2_vertices":r2_total,
            "n_triplet_neutral_uv_vertices":sum(
                v for k,v in counts.items() if k.startswith("neutral_UV")),
            "passes_basic_triplet_ct_gate":uv_total>0 and r2_total>0,
            "triplet_ct_vertex_breakdown":dict(counts),
            "status":"STRUCTURAL_PREFLIGHT_ONLY_NOT_PHYSICS_VALIDATED",
            "warnings":[
                "Static counterterm presence is not a check of correct UV renormalization.",
                "Require EW poles, FKS QED soft/collinear, gauge independence, SM limit.",
                "Require QED photon PDF (NNPDF40_nlo_as_01180_qed ID 335900).",
                "Keep mass splitting separate from hard-process tree parameters.",
            ]}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ufo",type=Path,
        default=Path("theory/dirac_triplet_ew/generated/"
                     "DiracTripletEW_QCDEW_CANDIDATE_UFO"))
    ap.add_argument("--json-out",type=Path)
    a=ap.parse_args()
    result=inspect_ufo(a.ufo)
    print(json.dumps(result,indent=2))
    if a.json_out:
        a.json_out.parent.mkdir(parents=True,exist_ok=True)
        a.json_out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    if not result["passes_basic_triplet_ct_gate"]:
        raise SystemExit("FAILED: BSM-specific UV+R2 are required for NLO EW.")
    print("PASS: BSM UV/R2 structural presence ONLY. Requires physical loop validation.")


if __name__=="__main__":
    main()
