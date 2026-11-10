"""Single-event CLI for independent quirk-string research prototype."""
import argparse
import json
from pathlib import Path
from .dynamics import Parameters, propagate


def main():
    p=argparse.ArgumentParser(description='QuirkSim straight-string pilot; no annihilation model')
    p.add_argument('--model',choices=['triplet','doublet'],default='triplet')
    p.add_argument('--mass',type=float,default=225)
    p.add_argument('--lambda-ev',type=float,default=774)
    p.add_argument('--beta',type=float,default=.4)
    p.add_argument('--impact-cm',type=float,default=0)
    p.add_argument('--duration-ps',type=float,default=500)
    p.add_argument('--seed',type=int,default=1)
    p.add_argument('--mode',choices=['stable','no_recoil','recoil'],default='recoil')
    p.add_argument('--trace',action='store_true',help='store all integrator trajectory points')
    p.add_argument('--max-step-ps',type=float,default=2.0)
    p.add_argument('--output',default='result.json')
    a=p.parse_args()
    params=Parameters(representation=a.model,mass_neutral_gev=a.mass,
        lambda_ev=a.lambda_ev,beta=a.beta,impact_cm=a.impact_cm,
        duration_ps=a.duration_ps,seed=a.seed,mode=a.mode,
        store_trace=a.trace,max_step_ps=a.max_step_ps)
    result=propagate(params)
    Path(a.output).parent.mkdir(parents=True,exist_ok=True)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('trace_time_ps_separation_cm','turning_points')},indent=2))


if __name__=='__main__':
    main()
