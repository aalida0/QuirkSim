"""Machine-readable exact radial NG benchmark, repeated-crossing conservation."""
import argparse
import json
from pathlib import Path
from quirksim.radial_ng_reference import RadialNGParameters, ExactRadialNG


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--duration-ps',type=float,default=1000.)
    p.add_argument('--mass',type=float,default=225.)
    p.add_argument('--lambda-ev',type=float,default=774.)
    p.add_argument('--initial-beta',type=float,default=.4)
    p.add_argument('--initial-separation-cm',type=float,default=.1)
    p.add_argument('--output',default='results/radial_ng_1000ps.json')
    a=p.parse_args()
    solver=ExactRadialNG(RadialNGParameters(
        mass_neutral_gev=a.mass,lambda_ev=a.lambda_ev,
        initial_beta=a.initial_beta,
        initial_separation_cm=a.initial_separation_cm))
    summary=solver.summary(a.duration_ps,samples=1001)
    summary['nominal_validation']='REDUCED_RADIAL_REFERENCE_ONLY'
    summary['explicit_crossing_times_ps']=[
        solver.first_crossing_ps+k*solver.crossing_interval_ps
        for k in range(summary['crossings_at_end'])
        if solver.first_crossing_ps+k*solver.crossing_interval_ps<=a.duration_ps]
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(summary,indent=2)+'\n')
    print('[radial-reference] '+json.dumps(summary),flush=True)


if __name__=='__main__':
    main()
