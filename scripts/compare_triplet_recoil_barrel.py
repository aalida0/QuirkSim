#!/usr/bin/env python3
"""Paired recoil-minus-no-recoil pixel-cylinder crossing diagnostics.

Compares EXACTLY matched production records and random seeds. The no-recoil
counterfactual is deliberately unphysical: mass and charge update on weak
decay but pion/lepton momentum is not subtracted from the endpoint.
Neither branch provides detector hit/reconstruction efficiencies.
"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path


def read(path):
    rows = {}
    with Path(path).open(encoding="utf-8") as f:
        for lineno, line in enumerate(f, 1):
            if not line.strip():
                continue
            record = json.loads(line)
            key = record.get("production_event_index")
            if key is None or key in rows:
                raise ValueError(f"Missing/duplicate production event ID in {path} at line {lineno}")
            rows[key] = record
    if not rows:
        raise ValueError(f"Empty results file: {path}")
    return rows


def first_charge_map(event):
    indexed = {}
    for row in event["first_outward_barrel_crossings"]:
        key = (int(row["endpoint"]), float(row["radius_cm"]))
        if key in indexed:
            raise ValueError(f"Duplicate outward crossing: {key}")
        indexed[key] = int(row["charge_at_crossing_e"]) != 0
    return indexed


def compare(recoil_path, no_recoil_path):
    reco = read(recoil_path)
    ctrl = read(no_recoil_path)
    if set(reco) != set(ctrl):
        raise ValueError("Unmatched production event IDs between runs")
    n = len(reco)
    radii = None
    event_pairs = {}
    endpoint_pairs = {}
    decay_pairs = Counter()

    for event_id in sorted(reco):
        a, b = reco[event_id], ctrl[event_id]
        if a.get("params", {}).get("mode") != "recoil" or b.get("params", {}).get("mode") != "no_recoil":
            raise ValueError(f"Wrong simulation modes for event {event_id}")
        for field in ("initial_pdg_ids", "initial_charges_e", "initial_p4_gev",
                      "initial_charged_endpoints", "pair_pt_gev"):
            if a.get(field) != b.get(field):
                raise ValueError(f"Unmatched production field {field} for event {event_id}")
        settings_a = dict(a["params"])
        settings_b = dict(b["params"])
        settings_a.pop("mode")
        settings_b.pop("mode")
        if settings_a != settings_b:
            raise ValueError(f"Runs differ in settings other than mode at event {event_id}")
        ra = tuple(map(float, a["barrel_layer_radii_cm"]))
        rb = tuple(map(float, b["barrel_layer_radii_cm"]))
        if ra != rb or (radii is not None and ra != radii):
            raise ValueError(f"Inconsistent layer geometry at event {event_id}")
        radii = ra
        crossed_a, crossed_b = first_charge_map(a), first_charge_map(b)
        decay_pairs[(a["n_weak_decays"], b["n_weak_decays"])] += 1
        for radius in radii:
            was_charged = [i for i,q in enumerate(a["initial_charges_e"]) if q!=0]
            # Endpoint paired events: each originally charged endpoint is one trial.
            for i in was_charged:
                key=(i,radius)
                cat=(bool(crossed_a.get(key,False)),bool(crossed_b.get(key,False)))
                endpoint_pairs.setdefault(radius,Counter())[cat]+=1
            a_event = any(crossed_a.get((i,radius),False) for i in was_charged)
            b_event = any(crossed_b.get((i,radius),False) for i in was_charged)
            event_pairs.setdefault(radius,Counter())[(a_event,b_event)] += 1

    def stats(counts, denom):
        both=counts[(True,True)]
        recoil_only=counts[(True,False)]
        control_only=counts[(False,True)]
        neither=counts[(False,False)]
        if both+recoil_only+control_only+neither!=denom:
            raise ValueError("Inconsistent paired category count")
        return {
            "denominator":denom,
            "both":both, "recoil_only":recoil_only,
            "no_recoil_only":control_only, "neither":neither,
            "recoil_fraction":(both+recoil_only)/denom,
            "no_recoil_fraction":(both+control_only)/denom,
            "paired_delta_fraction":(recoil_only-control_only)/denom,
        }

    results = []
    for radius in radii:
        ect=endpoint_pairs[radius]
        endpoint_denom=sum(ect.values())
        results.append({
            "radius_cm":radius,
            "event_level_at_least_one_charged_crossing":stats(event_pairs[radius],n),
            "endpoint_level_initially_charged_crossing":stats(ect,endpoint_denom),
        })
    return {
        "events":n, "radii":results,
        "decay_count_pairs":[{"recoil":k[0],"no_recoil":k[1],"events":count}
                             for k,count in sorted(decay_pairs.items())],
        "caveat":"No-recoil is an unphysical matched control; idealized infinite cylinders, no detector hits."
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("recoil",type=Path)
    p.add_argument("no_recoil",type=Path)
    p.add_argument("--json-out",type=Path)
    args=p.parse_args()
    out=compare(args.recoil,args.no_recoil)
    print(f"Matched production events: {out['events']}")
    print("R [cm]   recoil >=1  no_recoil >=1   paired delta    gained   lost")
    for row in out["radii"]:
        r=row["radius_cm"]
        q=row["event_level_at_least_one_charged_crossing"]
        print(f"{r:6.3f}     {q['recoil_fraction']:7.1%}     {q['no_recoil_fraction']:7.1%}"
              f"     {q['paired_delta_fraction']:+9.1%}      "
              f"{q['recoil_only']:4d}    {q['no_recoil_only']:4d}")
    print("Recoil-only = gains; no-recoil-only = losses.")
    print("No-recoil is unphysical; this estimates a controlled recoil sensitivity, NOT tracking efficiency.")
    if args.json_out:
        args.json_out.parent.mkdir(parents=True,exist_ok=True)
        args.json_out.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
        print("Wrote",args.json_out)


if __name__=="__main__":
    main()
