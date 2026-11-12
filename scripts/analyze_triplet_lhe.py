#!/usr/bin/env python3
"""Inspect the six LO hard-production channels in an unweighted MG5 LHE sample.

This reads the *hard LHE* event, not showered/decayed particles.
Each event must contain precisely two outgoing quirk endpoints.
Run from the QuirkSim checkout:
  python3 scripts/analyze_triplet_lhe.py results/pp_dirac_triplet_13p6TeV_LO/Events/run_01/unweighted_events.lhe.gz --sigma-pb 2.779 --n-ic 2
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import json
from math import atan2, hypot, log, sqrt
from pathlib import Path
from statistics import mean, median

PDGS = {9000111, 9000112, 9000113}
QCHARGE = {9000111: 1, -9000111: -1,
           9000112: 0, -9000112: 0,
           9000113: -1, -9000113: 1}
NAMES = {9000111: "q1+", -9000111: "q1-",
         9000112: "q0", -9000112: "q0~",
         9000113: "q2-", -9000113: "q2+"}
CHANNELS = {
    tuple(sorted((9000111, -9000111))): "charged_pair_q1",
    tuple(sorted((9000113, -9000113))): "charged_pair_q2",
    tuple(sorted((9000111, -9000112))): "associated_Wplus_q1",
    tuple(sorted((9000112, -9000113))): "associated_Wplus_q2",
    tuple(sorted((-9000111, 9000112))): "associated_Wminus_q1",
    tuple(sorted((-9000112, 9000113))): "associated_Wminus_q2",
}


def iter_lhe_endpoints(path: Path):
    """Yield a list of two (PDG,px,py,pz,E,m) tuples per LHE event."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        inside = False
        event_lines = []
        for line in stream:
            stripped = line.strip()
            if stripped.startswith("<event>"):
                if inside:
                    raise ValueError("Nested LHE <event> tags")
                inside, event_lines = True, []
            elif stripped.startswith("</event>"):
                if not inside:
                    raise ValueError("Unexpected </event>")
                payload = [s for s in event_lines if s and not s.startswith("#")]
                if not payload:
                    raise ValueError("Empty LHE event")
                try:
                    nup = int(payload[0].split()[0])
                except (ValueError, IndexError) as exc:
                    raise ValueError("Invalid LHE event header") from exc
                if len(payload) < 1 + nup:
                    raise ValueError("Truncated LHE event")
                endpoints = []
                for row in payload[1:1+nup]:
                    cols = row.split()
                    if len(cols) < 11:
                        raise ValueError("Incomplete LHE particle record: " + row)
                    pdg, status = int(cols[0]), int(cols[1])
                    if status == 1 and abs(pdg) in PDGS:
                        px, py, pz, energy, mass = map(float, cols[6:11])
                        endpoints.append((pdg, px, py, pz, energy, mass))
                if len(endpoints) != 2:
                    raise ValueError("Expected exactly two outgoing quirks, got %d" % len(endpoints))
                yield endpoints
                inside = False
            elif inside:
                event_lines.append(stripped)
        if inside:
            raise ValueError("Unterminated LHE event")


def analyze(path: Path, *, sigma_pb=None, n_ic=2):
    if n_ic < 1:
        raise ValueError("N_IC must be >= 1")
    if sigma_pb is not None and sigma_pb < 0:
        raise ValueError("Cross section must be nonnegative")
    counts = Counter()
    charged_counts = Counter()
    pt_charged, abseta_charged, pt_pair, invmass_pair = [], [], [], []
    n = 0
    for endpoints in iter_lhe_endpoints(path):
        n += 1
        key = tuple(sorted(p[0] for p in endpoints))
        name = CHANNELS.get(key)
        if name is None:
            raise ValueError("Unexpected quirk pair: " + repr([NAMES.get(p,p) for p in key]))
        counts[name] += 1
        charged = sum(QCHARGE[p[0]] != 0 for p in endpoints)
        charged_counts[charged] += 1
        px_sum = sum(p[1] for p in endpoints)
        py_sum = sum(p[2] for p in endpoints)
        pz_sum = sum(p[3] for p in endpoints)
        e_sum = sum(p[4] for p in endpoints)
        pt_pair.append(hypot(px_sum, py_sum))
        invmass_pair.append(sqrt(max(0.0, e_sum*e_sum - px_sum*px_sum
                                     - py_sum*py_sum - pz_sum*pz_sum)))
        for pdg, px, py, pz, energy, mass in endpoints:
            if QCHARGE[pdg] == 0:
                continue
            pt = hypot(px, py)
            pt_charged.append(pt)
            # Handle exactly forward endpoints without infinities.
            pabs = hypot(pt, pz)
            if pabs > abs(pz) and pt > 0:
                abseta_charged.append(abs(0.5*log((pabs+pz)/(pabs-pz))))
    if n == 0:
        raise ValueError("No LHE events found")
    channel_data = []
    for key, name in CHANNELS.items():
        k = counts[name]
        entry = {"name":name, "events":k, "fraction":k/n}
        if sigma_pb is not None:
            entry["sigma_pb_per_hidden_color"] = sigma_pb*k/n
            entry["sigma_pb_NIC_scaled"] = n_ic*sigma_pb*k/n
        channel_data.append(entry)
    return {
        "file":str(path), "events":n, "n_ic":n_ic,
        "sigma_pb_per_hidden_color":sigma_pb,
        "sigma_pb_NIC_scaled":None if sigma_pb is None else sigma_pb*n_ic,
        "channels":channel_data,
        "initial_charged_endpoints":{str(i):charged_counts[i] for i in (0,1,2)},
        "mean_charged_pt_GeV":mean(pt_charged) if pt_charged else None,
        "median_charged_pt_GeV":median(pt_charged) if pt_charged else None,
        "mean_abs_eta_charged":mean(abseta_charged) if abseta_charged else None,
        "mean_pair_pt_GeV":mean(pt_pair),
        "mean_pair_invariant_mass_GeV":mean(invmass_pair),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lhe", type=Path, help="MG5 unweighted_events.lhe[.gz]")
    parser.add_argument("--sigma-pb", type=float, default=None,
                        help="MG5 inclusive LO cross section (before N_IC scaling)")
    parser.add_argument("--n-ic", type=int, default=2)
    parser.add_argument("--json-out", type=Path, default=None)
    args = parser.parse_args()
    result = analyze(args.lhe, sigma_pb=args.sigma_pb, n_ic=args.n_ic)
    n = result["events"]
    print("LO Dirac-triplet hard-production LHE composition")
    print(f"Events: {n}; N_IC: {args.n_ic}")
    print(f"{'Channel':26s} {'N':>6s} {'Fraction':>11s}" +
          (f" {'sigma_LO [pb]':>16s}" if args.sigma_pb is not None else ""))
    for row in result["channels"]:
        label = f"{row['name']:26s} {row['events']:6d} {100*row['fraction']:10.2f}%"
        if args.sigma_pb is not None:
            label += f" {row['sigma_pb_NIC_scaled']:16.5f}"
        print(label)
    print("Initial charged endpoints:", result["initial_charged_endpoints"])
    for key in ("mean_charged_pt_GeV", "median_charged_pt_GeV",
                "mean_abs_eta_charged", "mean_pair_pt_GeV",
                "mean_pair_invariant_mass_GeV"):
        print(f"{key}: {result[key]:.4f}" if result[key] is not None else f"{key}: n/a")
    if args.sigma_pb is not None:
        print(f"Total LO sigma * N_IC = {result['sigma_pb_NIC_scaled']:.6g} pb")
        print("Channel cross sections estimated from fractions in an unweighted sample;")
        print("finite-sample fluctuations are NOT included in MG5 integration errors.")
    print("Note: 1 charged initial endpoint != 1 weak decay. No neutral-neutral DY.")
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
        print("Wrote", args.json_out)


if __name__ == "__main__":
    main()
