#!/usr/bin/env python3
"""Compare paired recoil/no-recoil post-decay closest turning points.

The reported observable is the smallest *post-first-decay* local minimum of
the separation within the integration window. A missing local minimum is
undefined/censored, never zero or infinite. This proxy is not a physical
annihilation probability or full string-worldsheet calculation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from math import sqrt, isfinite
from pathlib import Path
from statistics import mean, median


def read_records(path):
    out = {}
    with Path(path).open(encoding="utf-8") as f:
        for line_number, line in enumerate(f, 1):
            if not line.strip():
                continue
            record = json.loads(line)
            ident = record.get("production_event_index")
            if ident is None or ident in out:
                raise ValueError(f"Missing or duplicate event ID in {path}:{line_number}")
            out[ident] = record
    if not out:
        raise ValueError("Empty input: " + str(path))
    return out


def quantile(sorted_values, proportion):
    if not sorted_values:
        return None
    x = (len(sorted_values)-1)*proportion
    lo=int(x)
    hi=min(lo+1, len(sorted_values)-1)
    return sorted_values[lo]*(hi-x)+sorted_values[hi]*(x-lo)


def summary_metric(values):
    data = sorted(values)
    if not data:
        return {"n":0, "median":None, "mean":None,
                "q10":None, "q90":None}
    return {"n":len(data), "median":median(data), "mean":mean(data),
            "q10":quantile(data,.1), "q90":quantile(data,.9)}


def compare(recoil_file, no_recoil_file):
    arows, brows = read_records(recoil_file), read_records(no_recoil_file)
    if set(arows)!=set(brows):
        raise ValueError("Event IDs differ; cannot make paired comparison")
    categories=Counter()
    differences_um=[]
    relative_differences=[]
    momentum_differences=[]
    first_time_deltas_ps=[]
    decay_histories=Counter()
    details=[]
    for key in sorted(arows):
        a,b=arows[key],brows[key]
        if a.get("params",{}).get("mode")!="recoil" or b.get("params",{}).get("mode")!="no_recoil":
            raise ValueError(f"Wrong modes in event {key}")
        for field in ("initial_pdg_ids","initial_charges_e","initial_p4_gev",
                      "initial_charged_endpoints","pair_pt_gev"):
            if a.get(field)!=b.get(field):
                raise ValueError(f"Production mismatch {field} in event {key}")
        pa,pb=dict(a["params"]),dict(b["params"])
        pa.pop("mode")
        pb.pop("mode")
        if pa!=pb:
            raise ValueError(f"Unmatched simulation settings in event {key}")
        da,db=a["decays"],b["decays"]
        deca=[(d["endpoint"],d["channel"]) for d in da]
        decb=[(d["endpoint"],d["channel"]) for d in db]
        decay_histories[(len(da),len(db))]+=1
        if deca and decb and deca[0]!=decb[0]:
            raise ValueError(f"Different first weak decay in paired event {key}")
        if da and db:
            first_time_deltas_ps.append(float(da[0]["time_ps"])-float(db[0]["time_ps"]))
        ra=a.get("post_decay_closest_separation_cm")
        rb=b.get("post_decay_closest_separation_cm")
        if ra is not None and (not isfinite(float(ra)) or ra<0):
            raise ValueError(f"Invalid recoil closest approach in event {key}")
        if rb is not None and (not isfinite(float(rb)) or rb<0):
            raise ValueError(f"Invalid control closest approach in event {key}")
        if ra is None and rb is None:
            category="neither_defined"
        elif ra is None:
            category="control_only_defined"
        elif rb is None:
            category="recoil_only_defined"
        else:
            category="both_defined"
            delta_um=(float(ra)-float(rb))*1e4
            differences_um.append(delta_um)
            if float(rb)>0:
                relative_differences.append((float(ra)-float(rb))/float(rb))
        categories[category]+=1
        momenta_a=a.get("final_endpoint_p3_gev")
        momenta_b=b.get("final_endpoint_p3_gev")
        if momenta_a is not None and momenta_b is not None:
            if len(momenta_a)!=2 or len(momenta_b)!=2:
                raise ValueError(f"Invalid final momenta in event {key}")
            delta=sqrt(sum((float(momenta_a[i][j])-float(momenta_b[i][j]))**2
                           for i in range(2) for j in range(3)))
            momentum_differences.append(delta)
        details.append({
            "event_id":key,
            "initial_charged":a["initial_charged_endpoints"],
            "weak_decays_recoil":len(da),
            "weak_decays_control":len(db),
            "closest_recoil_um":None if ra is None else float(ra)*1e4,
            "closest_no_recoil_um":None if rb is None else float(rb)*1e4,
            "paired_difference_um":
                (float(ra)-float(rb))*1e4 if ra is not None and rb is not None else None,
            "turning_point_category":category,
        })
    return {
        "events":len(arows), "paired_post_decay_turning_point_availability":dict(categories),
        "paired_difference_recoil_minus_control_um":summary_metric(differences_um),
        "fraction_recoil_larger_conditional":
            sum(d>0 for d in differences_um)/len(differences_um) if differences_um else None,
        "fraction_recoil_smaller_conditional":
            sum(d<0 for d in differences_um)/len(differences_um) if differences_um else None,
        "fraction_nonzero_conditional":
            sum(d!=0 for d in differences_um)/len(differences_um) if differences_um else None,
        "paired_fractional_difference":summary_metric(relative_differences),
        "final_pair_momenta_euclidean_difference_gev":summary_metric(momentum_differences),
        "first_weak_decay_lab_time_recoil_minus_control_ps":summary_metric(first_time_deltas_ps),
        "decay_count_pairings": [
            {"recoil":a,"no_recoil":b,"events":v}
            for (a,b),v in sorted(decay_histories.items())
        ],
        "events_detail":details,
        "interpretation":{
            "positive_difference":"Recoil has a larger post-decay closest turning-point distance",
            "negative_difference":"Recoil has a smaller post-decay closest turning-point distance",
            "undefined":"No post-first-decay local turning point observed before t_max; never substitute 0",
            "scope":"Straight-string endpoint-force model; no worldsheet, radiation or annihilation",
            "no_recoil":"Unphysical reference used solely to isolate the recoil response",
        },
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("recoil",type=Path)
    p.add_argument("no_recoil",type=Path)
    p.add_argument("--json-out",type=Path)
    args=p.parse_args()
    report=compare(args.recoil,args.no_recoil)
    print("Matched events:",report["events"])
    print("Turning-point availability:",report["paired_post_decay_turning_point_availability"])
    print("Delta r_min = r_min(recoil) - r_min(no recoil) [um]:")
    print(json.dumps(report["paired_difference_recoil_minus_control_um"],indent=2))
    print("Fraction with larger closest separation (conditional):",
          report["fraction_recoil_larger_conditional"])
    print("Fraction with smaller closest separation (conditional):",
          report["fraction_recoil_smaller_conditional"])
    print("Delta first weak-decay lab time [ps]:",
          report["first_weak_decay_lab_time_recoil_minus_control_ps"])
    print("Final endpoint momentum-vector difference [GeV]:",
          report["final_pair_momenta_euclidean_difference_gev"])
    print("Caution: post-decay turning-point separation, not annihilation probability.")
    if args.json_out:
        args.json_out.parent.mkdir(parents=True,exist_ok=True)
        args.json_out.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
        print("Wrote",args.json_out)


if __name__=="__main__":
    main()
