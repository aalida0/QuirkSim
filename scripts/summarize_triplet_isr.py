#!/usr/bin/env python3
"""Summarize Pythia8 ISR-induced pair recoil and initial charge topologies."""
import argparse
from collections import Counter
import json
from pathlib import Path
import statistics

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('jsonl',type=Path)
    p.add_argument('--min-pt',type=float,default=20.)
    a=p.parse_args()
    counts=Counter()
    pts=[]
    n=0
    with a.jsonl.open() as f:
        for line in f:
            if not line.strip():continue
            rec=json.loads(line)
            if rec['stage']!='PYTHIA8_LO_QCD_ISR_FSR':
                raise ValueError('Expected showered Pythia8 event')
            if len(rec['endpoints'])!=2:
                raise ValueError('Need precisely 2 quirk endpoints')
            counts[sum(x['charge_e']!=0 for x in rec['endpoints'])]+=1
            pts.append(float(rec['pair_pt_gev']))
            n+=1
    if not n:raise ValueError('Empty Pythia8 sample')
    pts.sort()
    print('Pythia8 ordinary-QCD ISR/FSR (unmerged LO 2->2)')
    print('Events:',n)
    for key in (0,1,2):
        print(f'Initially charged endpoints {key}: {counts[key]} ({counts[key]/n:.1%})')
    print(f'Mean pair pT: {statistics.mean(pts):.3f} GeV')
    print(f'Median pair pT: {statistics.median(pts):.3f} GeV')
    print(f'90th percentile pair pT: {pts[max(0,int(.9*(n-1)))]:.3f} GeV')
    print(f'Fraction pair pT > {a.min_pt:g} GeV: {sum(x>a.min_pt for x in pts)/n:.1%}')
    print('Caution: no hidden radiation, hadronization, MPI, jet matching or detector.')
if __name__=='__main__':
    main()
