#!/usr/bin/env python3
"""Run QuirkSim on final quirk endpoint four-vectors from Pythia8 ISR JSONL."""
import argparse
import json
from pathlib import Path
from quirksim.production_dynamics import ProductionParameters, read_events, simulate_produced


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("jsonl",type=Path)
    p.add_argument("--nevents",type=int,default=3)
    p.add_argument("--duration-ps",type=float,default=200.)
    p.add_argument("--lambda-ev",type=float,default=774.)
    p.add_argument("--b-tesla",type=float,default=2.)
    p.add_argument("--mode",choices=['stable','no_recoil','recoil'],default='recoil')
    p.add_argument("--force-model",choices=['ng_endpoint','central_proxy'],default='ng_endpoint')
    p.add_argument("--seed",type=int,default=10)
    p.add_argument("--output",type=Path,default=Path('results/triplet_isr_propagated.jsonl'))
    a=p.parse_args()
    if a.nevents<=0:
        p.error('--nevents must be positive')
    cfg=ProductionParameters(lambda_ev=a.lambda_ev,b_tesla=a.b_tesla,
        duration_ps=a.duration_ps,mode=a.mode,force_model=a.force_model)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    counters={0:0,1:0,2:0}
    processed=0
    with a.output.open("w",encoding="utf-8") as dest:
        for i,event in enumerate(read_events(a.jsonl,limit=a.nevents)):
            this=ProductionParameters(**{**vars(cfg),'seed':a.seed+i})
            result=simulate_produced(event,this)
            counters[result['initial_charged_endpoints']]+=1
            dest.write(json.dumps(result)+"\n")
            processed+=1
            print(f"[QuirkSim ISR] {processed}/{a.nevents} "
                  f"initial_charge={result['initial_charged_endpoints']} "
                  f"decays={result['n_weak_decays']} "
                  f"pair_pt={result['pair_pt_gev']:.3f} GeV",flush=True)
    if processed==0:
        raise RuntimeError('No Pythia8 records found')
    print("Final:",json.dumps({'events':processed,'initial_charged_counts':counters,
                               'results':str(a.output)},indent=2))


if __name__=='__main__':
    main()
