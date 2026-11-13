#!/usr/bin/env python3
"""Streaming inspection of a potentially very large NLOCT QCD+EW .nlo file.

The full SM+triplet NLOCT file can exceed 100 MB. Do not call
Path.read_text() on the full file just to count UV/R2 entries.

Reports *structural* information (sections, approximate vertex record
counts, triplet field fingerprints and corrupt/incomplete output).
This cannot prove UV/IR cancellations or correct electroweak
renormalization. Do not substitute for FeynArts/MG5 physics checks.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

NLO_SECTIONS = {
    b"R2$vertlist": "R2",
    b"UV$vertlist": "UV",
    b"FR$CTparam": "CTparam",
}
TRIPLET_TOKENS = (
    b"Qp1", b"Qp1bar", b"Qzero", b"Qzerobar", b"Qm2", b"Qm2bar",
    b"F[61]", b"F[62]", b"F[63]",
)
DANGEROUS = (
    b"Indeterminate", b"ComplexInfinity", b"DirectedInfinity",
    b"$Failed", b"$Aborted",
)
ORDER_MARKER = b"FR$InteractionOrderPerturbativeExpansion"
CLASS_MARKER = b"CT$Model"
GEN_MARKER = b"CT$GenericModel"
MAX_SAMPLES = 4


def inspect(path: Path) -> dict:
    path=path.resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    section=None
    found=set()
    counts=Counter()
    maxline=Counter()
    largest={k:[] for k in NLO_SECTIONS.values()}
    triplet={k:Counter() for k in NLO_SECTIONS.values()}
    dangerous=Counter()
    metadata=Counter()
    complete=set()
    sha=hashlib.sha256()
    last_nonempty=b""
    total_lines=0
    nested_openers=Counter()
    for line in path.open("rb"):
        sha.update(line)
        total_lines+=1
        stripped=line.strip()
        if not stripped:
            continue
        last_nonempty=stripped
        for marker in (ORDER_MARKER,CLASS_MARKER,GEN_MARKER):
            if marker in line:
                metadata[marker.decode("ascii")]+=1
        for needle in DANGEROUS:
            if needle in line:
                dangerous[needle.decode("ascii")]+=line.count(needle)
        for marker,name in NLO_SECTIONS.items():
            if stripped.startswith(marker+b" = {") or re.match(re.escape(marker)+rb"\s*=\s*\{",stripped):
                if section is not None:
                    nested_openers[section]+=1
                found.add(name)
                section=name
                break
        if section is None:
            continue
        if any(stripped.startswith(marker) for marker in NLO_SECTIONS):
            continue
        if stripped==b"};" or stripped==b"}":
            complete.add(section)
            section=None
            continue
        counts[section]+=1
        maxline[section]=max(maxline[section],len(line))
        for token in TRIPLET_TOKENS:
            if token in line:
                triplet[section][token.decode("ascii")]+=line.count(token)
        if len(largest[section])<MAX_SAMPLES:
            largest[section].append({
                "line_number":total_lines,
                "line_length_bytes":len(line),
                "sha256":hashlib.sha256(line).hexdigest(),
                "prefix":line[:140].decode("utf-8","replace").replace("\n"," "),
            })
        # NLOCT writes the final list element as '<expr>};' on its line.
        if stripped.endswith(b"};"):
            complete.add(section)
            section=None
    missing=sorted(set(NLO_SECTIONS.values())-found)
    unclosed=sorted(set(NLO_SECTIONS.values())-complete)
    problems=[]
    if missing:
        problems.append("Required NLOCT sections absent: "+", ".join(missing))
    if unclosed:
        problems.append("Unclosed NLOCT section(s): "+", ".join(unclosed))
    if nested_openers:
        problems.append("Nested NLOCT list declarations: "+str(dict(nested_openers)))
    for name in ("UV","R2"):
        if counts[name]==0:
            problems.append(f"No vertex entry lines in {name} list")
    if not all(metadata.get(k.decode("ascii")) for k in (ORDER_MARKER,CLASS_MARKER,GEN_MARKER)):
        problems.append("Missing NLOCT scheme/FA provenance fields")
    if dangerous:
        problems.append("Non-finite or failed Wolfram token(s) in generated NLOCT file")
    if not last_nonempty.endswith(b"};"):
        problems.append("NLOCT file does not end with list closure '};'")
    has_triplet=any(triplet[k] for k in ("UV","R2"))
    if not has_triplet:
        # Do not flag as formal failure: FeynArts can renumber field
        # classes, and manual checks may need a mapping.
        problems.append("No recognizable triplet field names/class IDs in UV/R2; inspect class translation manually")
    return {
        "file":str(path),
        "bytes":path.stat().st_size,
        "sha256":sha.hexdigest(),
        "total_lines":total_lines,
        "sections_present":sorted(found),
        "sections_closed":sorted(complete),
        "vertex_record_line_counts_approximate":dict(counts),
        "longest_record_line_bytes":dict(maxline),
        "model_metadata_occurrences":dict(metadata),
        "triplet_field_occurrences":{k:dict(v) for k,v in triplet.items()},
        "representative_first_records":largest,
        "nonfinite_tokens":dict(dangerous),
        "findings":problems,
        "status":"CONTENT_STRUCTURAL_SCREEN_PASS_NOT_PHYSICS_CERTIFIED" if not problems
                 else "CONTENT_REQUIRES_REVIEW",
        "warning":"Structural screening cannot establish renormalization, UV/IR cancellation, gauge independence or N_IC=2 virtual-loop multiplicity.",
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nlo",type=Path,
        default=Path("theory/dirac_triplet_ew/generated/DiracTripletEW_QCDEW_CT.nlo"))
    ap.add_argument("--output",type=Path,default=Path("results/triplet_ew_nloct_content_report.json"))
    a=ap.parse_args()
    result=inspect(a.nlo)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    # Compact stdout: do not dump all representative vertex expressions.
    compact={key:result[key] for key in (
        "status","bytes","sha256","total_lines","sections_present",
        "sections_closed","vertex_record_line_counts_approximate",
        "longest_record_line_bytes","triplet_field_occurrences",
        "nonfinite_tokens","findings","warning")}
    print(json.dumps(compact,indent=2))
    print("Saved full structural report:",a.output)
    if result["status"]!="CONTENT_STRUCTURAL_SCREEN_PASS_NOT_PHYSICS_CERTIFIED":
        raise SystemExit("INSPECTION_REQUIRES_REVIEW: do not treat NLOCT output as certified.")


if __name__=="__main__":
    main()
