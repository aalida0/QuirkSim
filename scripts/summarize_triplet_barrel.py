#!/usr/bin/env python3
"""Summarize first outward crossings of infinite pixel-barrel cylinders.

The first crossing is NOT a detector hit: ignores finite barrel length,
sensor material, endcaps, charge deposition, reconstruction, and trigger.
The observation is right-censored at the simulation duration, usually 500 ps.
"""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
from statistics import mean,median


def summarize(path):
    total=0
    radius_data={}
    decay_radial=[]
    decay_times=[]
    with Path(path).open(encoding='utf-8') as infile:
        for line_no,line in enumerate(infile,1):
            if not line.strip():continue
            ev=json.loads(line)
            if 'first_outward_barrel_crossings' not in ev:
                raise ValueError('Old propagation record has no barrel geometry; rerun propagation')
            radii=ev.get('barrel_layer_radii_cm',[])
            if not radii:
                raise ValueError('Missing layer-radius specification')
            if not radius_data:
                radius_data={float(r):Counter() for r in radii}
            if set(radius_data)!={float(r) for r in radii}:
                raise ValueError('Inconsistent barrel layers across events')
            charges=ev['initial_charges_e']
            if len(charges)!=2:
                raise ValueError('Expected two endpoints')
            first={}
            for cross in ev['first_outward_barrel_crossings']:
                key=(int(cross['endpoint']),float(cross['radius_cm']))
                if key in first:
                    raise ValueError('Duplicate first crossing')
                first[key]=cross
            for r in radius_data:
                stats=radius_data[r]
                charged_in_event=0
                for i,original_charge in enumerate(charges):
                    crossing=first.get((i,r))
                    if original_charge!=0:
                        stats['initial_charged_endpoints']+=1
                        if crossing is None:
                            stats['charged_no_crossing_in_window']+=1
                        elif int(crossing['charge_at_crossing_e'])!=0:
                            stats['charged_at_first_crossing']+=1
                            charged_in_event+=1
                        else:
                            stats['neutral_at_first_crossing']+=1
                    else:
                        stats['initial_neutral_endpoints']+=1
                        if crossing is None:
                            stats['neutral_no_crossing_in_window']+=1
                        else:
                            if int(crossing['charge_at_crossing_e'])!=0:
                                raise ValueError('Initially neutral endpoint became charged')
                            stats['initial_neutral_crossed']+=1
                if charged_in_event>=1:stats['events_ge_one_charged_crossing']+=1
                if charged_in_event>=2:stats['events_ge_two_charged_crossings']+=1
            for decay in ev['decays']:
                rr=float(decay['radius_cm'])
                tt=float(decay['time_ps'])
                if not math.isfinite(rr) or rr<0 or not math.isfinite(tt):
                    raise ValueError('Invalid weak-decay vertex')
                decay_radial.append(rr)
                decay_times.append(tt)
            total+=1
    if total==0:
        raise ValueError('Empty propagated-event sample')
    rows=[]
    for r,stat in sorted(radius_data.items()):
        init=stat['initial_charged_endpoints']
        row={
            'radius_cm':r,'events':total,
            'initial_charged_endpoints':init,
            'charged_crossed_while_charged':stat['charged_at_first_crossing'],
            'charged_crossed_after_neutralizing':stat['neutral_at_first_crossing'],
            'initial_charged_not_crossed_in_window':stat['charged_no_crossing_in_window'],
            'fraction_initial_charged_crossed_while_charged':
                stat['charged_at_first_crossing']/init if init else None,
            'fraction_events_at_least_one_charged_crossing':
                stat['events_ge_one_charged_crossing']/total,
            'fraction_events_two_charged_crossings':
                stat['events_ge_two_charged_crossings']/total,
            'initial_neutral_endpoints':stat['initial_neutral_endpoints'],
            'initial_neutral_crossed':stat['initial_neutral_crossed'],
            'initial_neutral_not_crossed':stat['neutral_no_crossing_in_window'],
        }
        if (row['charged_crossed_while_charged']
                +row['charged_crossed_after_neutralizing']
                +row['initial_charged_not_crossed_in_window']!=init):
            raise ValueError('Charged-endpoint first-crossing categories do not sum')
        rows.append(row)
    return {
        'events':total,'layers':rows,
        'weak_decays_observed':len(decay_radial),
        'mean_decay_radius_cm':mean(decay_radial) if decay_radial else None,
        'median_decay_radius_cm':median(decay_radial) if decay_radial else None,
        'mean_decay_time_ps':mean(decay_times) if decay_times else None,
        'caveat':'Infinite barrel cylinders, no finite-z/endcap acceptance or hit simulation; finite propagation window.'
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('jsonl',type=Path)
    p.add_argument('--json-out',type=Path)
    a=p.parse_args()
    report=summarize(a.jsonl)
    print('Events:',report['events'],'Observed weak decays:',report['weak_decays_observed'])
    print('Mean decay radius [cm]:',report['mean_decay_radius_cm'])
    print('Median decay radius [cm]:',report['median_decay_radius_cm'])
    print('Radius [cm] | charged crossing / initially charged | after-neutral crossing | not reached | events >=1 charged cross (%)')
    for row in report['layers']:
        print(f"{row['radius_cm']:8.3f} | "
              f"{row['charged_crossed_while_charged']:4d}/{row['initial_charged_endpoints']:4d} "
              f"({100*row['fraction_initial_charged_crossed_while_charged']:5.1f}%) | "
              f"{row['charged_crossed_after_neutralizing']:4d} | "
              f"{row['initial_charged_not_crossed_in_window']:4d} | "
              f"{100*row['fraction_events_at_least_one_charged_crossing']:5.1f}%")
    print('CAUTION:',report['caveat'])
    if a.json_out:
        a.json_out.parent.mkdir(parents=True,exist_ok=True)
        a.json_out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('Wrote',a.json_out)
if __name__=='__main__':
    main()
