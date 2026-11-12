#!/usr/bin/env python3
"""Extract the standard 7-point NLO-QCD scale envelope from MG5_aMC 3.5.1.

This reuses the completed run's weighted scale_pdf_dependence.dat tables.
MG5 puts renormalization-factor indices fastest, factorization indices
slowest, and ensures the (muR,muF)=(1,1) pair is the first table entry.
The 7-point prescription excludes (0.5,2.0) and (2.0,0.5).

Crucial: use this BEFORE launching LO or another integration in the same
MG5 process directory, which may overwrite SubProcesses/job_status.pkl.

Read only trusted local MadGraph pickle files (pickle is not safe for
untrusted input).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import pickle
import re

from extract_triplet_nlo_pdf_replicas import (
    pdf_replica_summary,
    read_central_from_log,
    resolve_job_directory,
)

FACTORS = frozenset({0.5, 1.0, 2.0})


def get_scale_factors(run_card: Path, key: str) -> list[float]:
    """Read original run card and mimic MG5's central-factor normalization."""
    source = run_card.read_text(encoding="utf-8")
    m = re.search(
        r"^\s*([^\n=!]+?)\s*=\s*" + re.escape(key) + r"\b",
        source, re.MULTILINE | re.IGNORECASE,
    )
    if not m:
        raise ValueError(f"Missing '{key}' in {run_card}")
    vals = [float(v) for v in re.findall(
        r"(?<![A-Za-z0-9_.])[+-]?(?:\d+\.\d*|\.\d+|\d+)"
        r"(?:[eEdD][+-]?\d+)?",
        m.group(1).replace("D", "E").replace("d", "e"),
    )]
    if len(vals) != 3 or set(vals) != FACTORS:
        raise ValueError(f"Expect exactly 0.5,1,2 in {key}, found {vals}")
    # MG5 3.5.1 RunCard.check_validity(): the 1.0 factor is swapped first.
    if vals[0] != 1.0:
        idx = vals.index(1.0)
        vals[0], vals[idx] = vals[idx], vals[0]
    return vals


def read_scale_block(path: Path, expected_n: int = 9) -> tuple[int, list[float]]:
    """Read first 'scale variations' block, one metadata and one values line.

    The Fortran writer iterates outer muF and inner muR, using the
    (possibly normalized) RunCard rw_fscale and rw_rscale vectors.
    """
    with path.open(encoding="utf-8") as f:
        for line in f:
            if "scale variations:" in line.lower():
                head, raw = next(f, None), next(f, None)
                if head is None or raw is None:
                    raise ValueError(f"Incomplete scale block: {path}")
                header = head.split()
                if len(header) < 3:
                    raise ValueError(f"Invalid scale metadata: {head!r}")
                try:
                    dynamic, nR, nF = map(int, header[:3])
                    vals = [float(s.replace("D", "E").replace("d", "e"))
                            for s in raw.split()]
                except ValueError as exc:
                    raise ValueError(f"Bad scale numeric data in {path}") from exc
                if nR * nF != expected_n or len(vals) != expected_n:
                    raise ValueError(f"Expected {expected_n} scale weights in {path}; "
                                     f"header nR={nR}, nF={nF}, values={len(vals)}")
                if any(not math.isfinite(v) for v in vals):
                    raise ValueError(f"Non-finite scale weight in {path}")
                return dynamic, vals
    raise ValueError(f"Missing 'scale variations:' block: {path}")


def scale_envelope(
    process: Path,
    expected_central: float | None = None,
    relative_tolerance: float = 0.01,
) -> dict:
    process = process.resolve()
    # MG5 reweighting files reflect the most recent integration, so
    # use its current run card, not a historic banner or generated card.
    card = process / "Cards" / "run_card.dat"
    rvals = get_scale_factors(card, "rw_rscale")
    fvals = get_scale_factors(card, "rw_fscale")
    with (process / "SubProcesses" / "job_status.pkl").open("rb") as stream:
        jobs = pickle.load(stream)
    if not isinstance(jobs, list) or not jobs:
        raise ValueError("Missing MG5 integration jobs")
    totals = [0.] * 9
    dynamic_choice = None
    for j in jobs:
        if "wgt_frac" not in j:
            raise ValueError("Missing job weight fraction")
        wgt_frac = float(j["wgt_frac"])
        if not math.isfinite(wgt_frac):
            raise ValueError("Non-finite job weight fraction")
        path = resolve_job_directory(process, j) / "scale_pdf_dependence.dat"
        dyn, vals = read_scale_block(path)
        if dynamic_choice is None:
            dynamic_choice = dyn
        if dynamic_choice != dyn:
            raise ValueError("Inconsistent dynamic-scale choices across jobs")
        for i, val in enumerate(vals):
            totals[i] += wgt_frac * val
    if dynamic_choice != -1:
        raise ValueError(f"Expected dynamical_scale_choice=-1, got {dynamic_choice}")

    grid = {
        (rvals[r], fvals[f]): totals[f*len(rvals)+r]
        for f in range(len(fvals))
        for r in range(len(rvals))
    }
    if len(grid) != 9:
        raise ValueError("Nonunique scale factors")
    central = grid[(1.,1.)]
    if central <= 0 or (expected_central is not None and
                       abs(central-expected_central)/expected_central > relative_tolerance):
        raise ValueError(
            f"Scale center {central:.8g} pb mismatches expected "
            f"{expected_central} pb; check job_status and run_card. No uncertainty quoted."
        )
    seven = {
        pair: sig for pair,sig in grid.items()
        if pair not in {(0.5,2.0),(2.0,0.5)}
    }
    if len(seven) != 7:
        raise ValueError("Incorrect 7-point scale selection")
    low = min(seven.values())
    high = max(seven.values())
    return {
        "center_pb": central,
        "seven_low_pb": low,
        "seven_high_pb": high,
        "seven_down_pb": central-low,
        "seven_up_pb": high-central,
        "seven_down_percent": (central-low)/central*100,
        "seven_up_percent": (high-central)/central*100,
        "nine_down_percent": (central-min(grid.values()))/central*100,
        "nine_up_percent": (max(grid.values())-central)/central*100,
        "points_pb": {f"muR={p[0]:g},muF={p[1]:g}": v
                       for p,v in sorted(grid.items())},
        "njobs": len(jobs),
        "dynamic": dynamic_choice,
    }


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--process",type=Path,
        default=Path("results/pp_dirac_triplet_13p6TeV_NLO_NNPDF40"))
    parser.add_argument("--run-log",type=Path,
        default=Path("results/triplet_nlo_uncertainties.log"))
    parser.add_argument("--nic",type=int,default=2)
    parser.add_argument("--json-out",type=Path)
    a=parser.parse_args()
    if a.nic < 1:
        parser.error("--nic must be positive")
    expected=read_central_from_log(a.run_log)
    result=scale_envelope(a.process,expected_central=expected)
    pdf=pdf_replica_summary(a.process,expected_central=expected)
    if abs(result["center_pb"]-pdf["central_pb"])/pdf["central_pb"] > 0.01:
        raise ValueError("PDF and scale reweighting central values inconsistent")
    print(f"MG5 jobs: {result['njobs']}, dynamical scale: {result['dynamic']}")
    print(f"Central sigma (scale weights): {result['center_pb']:.7f} pb")
    print("Scale variation grid (muR/mu0,muF/mu0):")
    for p,value in result["points_pb"].items():
        print(f"  {p}: {value:.7f} pb")
    print(f"9-point envelope (MG5 cross-check): "
          f"+{result['nine_up_percent']:.3f}%/-{result['nine_down_percent']:.3f}%")
    print(f"7-point envelope: +{result['seven_up_percent']:.3f}%/"
          f"-{result['seven_down_percent']:.3f}%")
    print(f"7-point absolute: +{result['seven_up_pb']:.6f}/"
          f"-{result['seven_down_pb']:.6f} pb")
    print(f"For N_IC={a.nic}: {result['center_pb']*a.nic:.7f} pb "
          f"(7-point +{result['seven_up_pb']*a.nic:.6f}/"
          f"-{result['seven_down_pb']*a.nic:.6f} pb)")
    print(f"PDF-only (for N_IC={a.nic}): "
          f"+/-{pdf['pdf_sigma_pb']*a.nic:.6f} pb "
          f"({pdf['pdf_sigma_percent']:.3f}%)")
    print("Note: scale and PDF errors are separate; no LO K-factor computed.")
    if a.json_out is not None:
        a.json_out.parent.mkdir(parents=True,exist_ok=True)
        payload={"process":str(a.process.resolve()),"scale":result,
                 "pdf":pdf,"nic":a.nic,"source":"MG5 3.5.1 saved job weights"}
        a.json_out.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
        print("Saved:",a.json_out)


if __name__=="__main__":
    main()
