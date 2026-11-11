"""Aggregate near-contact shape refinement; distinguish tolerance vs mesh errors.

This script never labels a few-grid comparison as certified continuum NG.
The core-scale threshold is a diagnostic reference only, NOT annihilation.
"""
import argparse
import json
from pathlib import Path
from quirksim.dynamics import HBARC_GEV_CM


def audit(root, output):
    files=sorted(Path(root).glob('near_contact_*.json'))
    rows=[json.loads(p.read_text()) for p in files]
    reference=HBARC_GEV_CM/(774e-9)
    failed=[x for x in rows if x.get('status')!='PHYSICAL_ARCLENGTH_CUTOFF']
    by_mesh={}
    for r in rows:
        if r.get('status')=='PHYSICAL_ARCLENGTH_CUTOFF':
            by_mesh.setdefault(r['nodes'],[]).append(r)
    ordered={n:sorted(values,key=lambda r:r['numerical_max_step_ps'],reverse=True)
             for n,values in by_mesh.items()}
    base={}
    for n,rr in ordered.items():
        # Default h=.5ps, otherwise widest step.
        base[n]=min(rr,key=lambda r:abs(r['numerical_max_step_ps']-.5))
    mesh_sizes=sorted(base)
    compare=[]
    for na,nb in zip(mesh_sizes[:-1],mesh_sizes[1:]):
        a,b=base[na],base[nb]
        compare.append({
            'nodes':(na,nb),
            'chord_difference_cm':abs(a['endpoint_separation_cm']-b['endpoint_separation_cm']),
            'chord_difference_in_infracolor_lengths':
                abs(a['endpoint_separation_cm']-b['endpoint_separation_cm'])/reference,
            'bend_angle_difference_rad':abs(
                a['maximum_adjacent_segment_bend_angle_rad']-
                b['maximum_adjacent_segment_bend_angle_rad'])
        })
    time_step=[]
    for n,rr in ordered.items():
        if len(rr)<2:
            continue
        top=max(rr,key=lambda r:r['numerical_max_step_ps'])
        low=min(rr,key=lambda r:r['numerical_max_step_ps'])
        time_step.append({
           'nodes':n,'coarse_step_ps':top['numerical_max_step_ps'],
           'fine_step_ps':low['numerical_max_step_ps'],
           'closest_chord_difference_cm':abs(
              top['endpoint_separation_cm']-low['endpoint_separation_cm']),
           'difference_in_infracolor_lengths':
              abs(top['endpoint_separation_cm']-low['endpoint_separation_cm'])/reference
        })
    # Even if the last two resolutions agree, we still require a third
    # independent resolution to show that their differences are shrinking.
    last_three=mesh_sizes[-3:]
    demonstrated=False
    if len(last_three)==3:
        diffs=[abs(base[last_three[i]]['endpoint_separation_cm']-
                   base[last_three[i+1]]['endpoint_separation_cm']) for i in range(2)]
        demonstrated=(diffs[1]<diffs[0] and diffs[1]<reference*.1 and
                      all(x['difference_in_infracolor_lengths']<.1
                          for x in time_step))
    report={
       'status':('PREPINCH_FINE_GRID_CHORD_DIAGNOSTIC_PASS_NOT_CONTINUUM'
                 if demonstrated and not failed else
                 'BLOCKED_INFRACOLOR_SCALE_CHORD_CONVERGENCE'),
       'infracolor_reference_length_cm':reference,
       'file_count':len(rows),
       'mesh_sizes':mesh_sizes,
       'failures':failed,
       'mesh_comparisons':compare,
       'time_step_comparisons':time_step,
       'finest_resolution':base[mesh_sizes[-1]] if mesh_sizes else None,
       'limits':[
          'no NG worldsheet continuation at pinch',
          'single 2T pre-stretched slice, not pp production kinematics',
          'weak decays not coupled to new zero-mass gauge',
          'underconstrained string short-distance dynamics; even PASS is diagnostic',
          'geometric 1/Lambda is NOT annihilation contact radius']
    }
    out=Path(output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+'\n')
    print('[near-contact-audit] '+json.dumps(report),flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--input-dir',default='results')
    parser.add_argument('--output',default='results/near_contact_audit.json')
    args=parser.parse_args()
    audit(args.input_dir,args.output)
