#!/usr/bin/env python3
"""Targeted, memory-bounded classification of NLOCT UV and R2 external legs.

NLOCT writes one pair {external_fields, expression} per line. This
script reads each line separately (large coefficients can be ~2.5 MB),
extracts the FIRST top-level item (external fields), and never prints
the large amplitudes. Looking for triplet symbols anywhere in an entire
counterterm line is NOT sufficient: a triplet mass may occur inside an
SM-only coefficient. An explicit triplet external-leg audit is needed.

This is a data-forensics diagnostic and not proof of UV/IR cancellation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

SECTION_TOKENS={b"R2$vertlist": "R2",b"UV$vertlist": "UV"}
NEEDLES=(b"Qp1",b"Qzero",b"Qm2",b"Qp1bar",b"Qzerobar",b"Qm2bar")
MAX_DISTINCT=10000
SAMPLE_LIMIT=12


def extract_first_item(data: bytes) -> bytes | None:
    """Return the first subexpression in {first, rest} without parsing rest.

    Supports nested brackets/parentheses, strings and Mathematica's
    functional syntax. Does not interpret any input as Mathematica code.
    """
    s=data.strip()
    if not s.startswith(b"{"):
        return None
    stack=[]
    matching={ord("}"):ord("{"),ord("]"):ord("["),ord(")"):ord("(")}
    quoted=False
    escaped=False
    for i,c in enumerate(s):
        if quoted:
            if escaped:
                escaped=False
            elif c==ord("\\"):
                escaped=True
            elif c==ord('"'):
                quoted=False
            continue
        if c==ord('"'):
            quoted=True
            continue
        if c in (ord("{"),ord("["),ord("(")):
            stack.append(c)
            continue
        if c in matching:
            if not stack or stack[-1]!=matching[c]:
                return None
            stack.pop()
            if not stack:
                return None
            continue
        if c==ord(",") and len(stack)==1:
            return s[1:i].strip()
    return None


def inspect(path: Path, *, samples: int=SAMPLE_LIMIT) -> dict:
    path=path.resolve()
    section=None
    counts=Counter()
    header_counts={section:Counter() for section in SECTION_TOKENS.values()}
    explicit_external={section:Counter() for section in SECTION_TOKENS.values()}
    strange=[]
    examples={section:[] for section in SECTION_TOKENS.values()}
    sha=hashlib.sha256()
    for line_number,line in enumerate(path.open("rb"),1):
        sha.update(line)
        clean=line.strip()
        for key,name in SECTION_TOKENS.items():
            if clean.startswith(key+b" = {"):
                section=name
                break
        else:
            if section is None:
                continue
            if clean==b"};":
                section=None
                continue
            header=extract_first_item(clean)
            if header is None:
                if len(strange)<12:
                    strange.append({"line":line_number,"section":section,
                                    "prefix":clean[:160].decode("utf-8","replace")})
                counts[f"{section}_unparsed"]+=1
                if clean.endswith(b"};"):
                    section=None
                continue
            counts[section]+=1
            key=header.decode("utf-8","replace").replace("\n"," ")
            if len(header_counts[section])<MAX_DISTINCT or key in header_counts[section]:
                header_counts[section][key]+=1
            else:
                counts[f"{section}_distinct_header_limit_exceeded"]+=1
            for needle in NEEDLES:
                if needle in header:
                    explicit_external[section][needle.decode()]+=1
            if len(examples[section])<samples:
                examples[section].append({"line":line_number,
                                          "external_legs":key[:300]})
            if clean.endswith(b"};"):
                section=None
    nloct={"input":str(path),"file_sha256":sha.hexdigest(),
          "records_parsed":dict(counts),
          "explicit_triplet_external_header_counts":
              {k:dict(v) for k,v in explicit_external.items()},
          "distinct_external_legs_counts":
              {k:len(v) for k,v in header_counts.items()},
          "most_frequent_headers":
              {k:[{"header":h[:300],"count":n} for h,n in v.most_common(25)]
               for k,v in header_counts.items()},
          "initial_header_examples":examples,
          "unparsed_examples":strange,
          "warning":"Header tokens may be renumbered by FeynArts. A clean result is NOT a physics certification."}
    nloct["status"]=(
       "HEADERS_PARSED_REQUIRE_PHYSICS_VALIDATION" if
       counts.get("UV",0)>0 and counts.get("R2",0)>0 and
       counts.get("UV_unparsed",0)==0 and counts.get("R2_unparsed",0)==0
       else "HEADERS_NEED_MANUAL_REVIEW")
    return nloct


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nlo",type=Path,
        default=Path("theory/dirac_triplet_ew/generated/DiracTripletEW_QCDEW_CT.nlo"))
    p.add_argument("--output",type=Path,
        default=Path("results/triplet_ew_nloct_header_audit.json"))
    args=p.parse_args()
    report=inspect(args.nlo)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:report[k] for k in (
        "status","file_sha256","records_parsed",
        "explicit_triplet_external_header_counts",
        "distinct_external_legs_counts","most_frequent_headers",
        "initial_header_examples","unparsed_examples","warning")},indent=2))
    print("Report written:",args.output)
    if report["status"]!="HEADERS_PARSED_REQUIRE_PHYSICS_VALIDATION":
        raise SystemExit("NLOCT external-leg headers need manual review; no UFO export is warranted.")


if __name__=="__main__":
    main()
