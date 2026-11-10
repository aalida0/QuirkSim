"""Audit near-contact shape convergence without crossing the NG pinch.

One invocation = one resolution/step-size numerical experiment.  Results
are machine-readable and explicitly BLOCKED if a mesh-quality or KKT
condition prevents reaching the common geometric arclength threshold.
"""
import argparse
import json
from pathlib import Path
from equidistributed_gauge_stress import case
from quirksim.dynamics import HBARC_GEV_CM


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--nodes',type=int,required=True)
    p.add_argument('--max-step-ps',type=float,default=.5)
    p.add_argument('--rtol',type=float,default=2e-9)
    p.add_argument('--atol',type=float,default=1e-10)
    p.add_argument('--pinch-fraction',type=float,default=2e-5)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    if not 0<a.pinch_fraction<1 or a.max_step_ps<=0:
        p.error('invalid time step or arclength cutoff')
    row=case(a.nodes,2.,.25,1000.,a.pinch_fraction,
             max_step_ps=a.max_step_ps,rtol=a.rtol,atol=a.atol)
    ref=HBARC_GEV_CM/(774e-9)
    row['infracolor_reference_length_cm']=ref
    if row.get('status')=='PHYSICAL_ARCLENGTH_CUTOFF':
        row['chord_deficit_in_reference_lengths']=row['chord_deficit_cm']/ref
    row['physics_gate']='PRE_PINCH_ONLY_NO_CONTINUATION'
    out=Path(a.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(row,indent=2)+'\n')
    print('[near-contact-shape] '+json.dumps(row),flush=True)


if __name__=='__main__':
    main()
