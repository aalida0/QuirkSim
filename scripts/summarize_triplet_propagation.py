#!/usr/bin/env python3
"""Summarize QuirkSim initial-charge -> weak-decay -> final-charge transitions.

Reads scripts/propagate_triplet_isr.py output JSONL, with one record per
original Pythia8 event. Does NOT infer detection or MET efficiencies.
"""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
from statistics import mean, median


def summarize(path):
    matrix = Counter()
    by_initial = Counter()
    by_final = Counter()
    by_n_decay = Counter()
    transitions = defaultdict(list)
    decay_channels = Counter()
    decay_times_ps = []
    pt_values = []
    n = 0
    for line_no, line in enumerate(Path(path).open(encoding='utf-8'), 1):
        if not line.strip():
            continue
        rec = json.loads(line)
        initial = int(rec['initial_charged_endpoints'])
        final = int(rec['final_charged_endpoints'])
        decays = rec['decays']
        nd = int(rec['n_weak_decays'])
        if initial not in (1, 2) or final < 0 or final > initial:
            raise ValueError(f'Invalid initial/final topology at line {line_no}')
        if nd != len(decays) or initial - nd != final:
            raise ValueError(f'Weak-decay/charge-count mismatch at line {line_no}')
        pt = float(rec['pair_pt_gev'])
        if not math.isfinite(pt) or pt < 0:
            raise ValueError(f'Invalid pair pT at line {line_no}')
        for d in decays:
            if d['channel'] not in ('pion', 'electron_neutrino', 'muon_neutrino'):
                raise ValueError(f'Invalid weak-decay channel at line {line_no}')
            time = float(d['time_ps'])
            if not (0 <= time <= float(rec['end_time_ps'])+1e-6):
                raise ValueError(f'Invalid decay time at line {line_no}')
            decay_times_ps.append(time)
            decay_channels[d['channel']] += 1
        n += 1
        key = (initial, nd, final)
        matrix[key] += 1
        transitions[initial].append(nd)
        by_initial[initial] += 1
        by_final[final] += 1
        by_n_decay[nd] += 1
        pt_values.append(pt)
    if not n:
        raise ValueError('No propagated JSONL events found')
    rows = []
    for key, count in sorted(matrix.items()):
        ini, decays, fin = key
        rows.append({
            'initial_charged': ini, 'n_weak_decays': decays,
            'final_charged': fin, 'events': count,
            'fraction_all': count/n,
            'fraction_given_initial': count/by_initial[ini],
        })
    return {
        'events': n,
        'transitions': rows,
        'initial_charged_counts': {str(i): by_initial[i] for i in (1, 2)},
        'weak_decay_counts': {str(i): by_n_decay[i] for i in (0, 1, 2)},
        'final_charged_counts': {str(i): by_final[i] for i in (0, 1, 2)},
        'decay_channels': dict(decay_channels),
        'mean_decay_time_ps': mean(decay_times_ps) if decay_times_ps else None,
        'median_decay_time_ps': median(decay_times_ps) if decay_times_ps else None,
        'mean_pair_pt_gev': mean(pt_values),
        'median_pair_pt_gev': median(pt_values),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('jsonl', type=Path)
    parser.add_argument('--json-out', type=Path)
    args = parser.parse_args()
    report = summarize(args.jsonl)
    print(f"Events analyzed: {report['events']}")
    print(f"{'Initial':>8} {'Decays':>8} {'Final':>8} {'Events':>8} {'All (%)':>10} {'Within initial (%)':>20}")
    for row in report['transitions']:
        print(f"{row['initial_charged']:>8d} {row['n_weak_decays']:>8d} "
              f"{row['final_charged']:>8d} {row['events']:>8d} "
              f"{100*row['fraction_all']:>10.1f} {100*row['fraction_given_initial']:>20.1f}")
    print('Initial charged counts:', report['initial_charged_counts'])
    print('Number of weak decays:', report['weak_decay_counts'])
    print('Final charged counts:', report['final_charged_counts'])
    print('Weak decay channels:', report['decay_channels'])
    for key in ('mean_pair_pt_gev','median_pair_pt_gev','mean_decay_time_ps','median_decay_time_ps'):
        print(key, report[key])
    print('Caution: no track acceptance, annihilation, hadronization or full string worldsheet.')
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        print('Wrote:', args.json_out)


if __name__ == '__main__':
    main()
