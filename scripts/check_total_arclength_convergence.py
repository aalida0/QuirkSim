"""Fail-closed 2T triplet convergence at a COMMON TOTAL string arclength.

This is an independent CLI validation gate over N=5/7/11; it must NOT
mistake successful GitHub job execution for physical 1000ps completion.
All cases evolve until the earliest near-zero-length boundary, not through
an unvalidated folded/reconnected/annihilating worldsheet.
"""
import argparse
import json
import math
from pathlib import Path
from equidistributed_gauge_stress import case
from quirksim.dynamics import HBARC_GEV_CM


def validate(cutoff_fractions=(1e-3,4e-4), horizon_ps=1000,
             output_path='results/constrained_ng_common_arclength.json'):
    trials=[]
    for fraction in cutoff_fractions:
        results=[]
        for n in (5,7,11):
            row=case(n,2.,.25,horizon_ps,pinch_fraction=fraction)
            results.append(row)
            print('[common-arclength-case] '+json.dumps(row),flush=True)
        times=[r.get('time_ps',float('nan')) for r in results]
        residuals=[r.get('absolute_total_arclength_cutoff_residual_cm',
                          float('nan')) for r in results]
        energies=[abs(r.get('energy_drift_gev',float('nan'))) for r in results]
        chords=[r.get('endpoint_separation_cm',float('nan')) for r in results]
        # Geometric core reference only, NOT an annihilation radius.
        reference_core_cm=HBARC_GEV_CM/(774.0e-9)
        chord_spread=max(chords)-min(chords)
        constraints=[abs(r.get('max_equal_length_constraint_cm',float('nan')))
                     for r in results]
        passed=(
            all(r.get('status')=='PHYSICAL_ARCLENGTH_CUTOFF' and
                r.get('stops')==['total_arclength_cutoff'] for r in results)
            and all(math.isfinite(t) for t in times)
            and max(times)-min(times)<5e-4
            and all(math.isfinite(d) and d<2e-8 for d in residuals)
            and all(math.isfinite(e) and e<1e-6 for e in energies)
            and all(math.isfinite(g) and g<1e-6 for g in constraints)
            and all(r.get('gauge_point_mass_gev')==0 for r in results)
        )
        record={
            'cutoff_fraction':fraction,
            'common_total_cutoff_cm':0.1*fraction,
            'max_stop_time_spread_ps':max(times)-min(times),
            'max_arclength_event_residual_cm':max(residuals),
            'max_energy_drift_gev':max(energies),
            'max_endpoint_chord_spread_cm':chord_spread,
            'reference_infracolor_length_cm':reference_core_cm,
            'chord_spread_in_units_of_1_over_lambda':chord_spread/reference_core_cm,
            'core_scale_geometry_gate':(
                 'PASS_10PCT_REFERENCE_ONLY'
                 if chord_spread<.1*reference_core_cm else
                 'BLOCKED_INSUFFICIENT_CORE_SCALE_CONVERGENCE'),
            'max_equal_arclength_constraint_cm':max(constraints),
            'convergence_gate_pass':bool(passed),
            'results':results
        }
        trials.append(record)
        print('[common-arclength-check] '+json.dumps(
            {k:v for k,v in record.items() if k!='results'}),flush=True)
    report={'status':'PRE_PINCH_CONVERGENCE_VALIDATED_ONLY'
                      if all(t['convergence_gate_pass'] for t in trials)
                      else 'BLOCKED_MESH_REFINEMENT',
            'horizon_requested_ps':horizon_ps,
            'trials':trials,
            'warnings':[
              'STOP at near-zero arclength; 1000 ps has not been completed',
              'engineering cutoff, not annihilation or physical core radius',
              'same-point collider production not implemented',
              'no weak recoils in this mesh stress test',
              'no continuation of dynamically moving string through contact',
              'nonzero solenoid: mechanical momentum changes under magnetic force'
            ]}
    path=Path(output_path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2)+'\n')
    print('[common-arclength-summary] '+json.dumps({
       'status':report['status'],
       'trials':[{'cutoff':t['common_total_cutoff_cm'],
                  'spread_ps':t['max_stop_time_spread_ps'],
                  'passed':t['convergence_gate_pass'],
                  'core_geometry_gate':t['core_scale_geometry_gate']}
                  for t in trials]}),flush=True)
    if report['status']!='PRE_PINCH_CONVERGENCE_VALIDATED_ONLY':
        raise AssertionError('common-arclength mesh convergence validation failed')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--fractions',default='0.001,0.0004',
                        help='comma-separated total-arclength/initial-length thresholds')
    parser.add_argument('--output',default='results/constrained_ng_common_arclength.json')
    args=parser.parse_args()
    fractions=tuple(float(z.strip()) for z in args.fractions.split(','))
    if not fractions or any(not 0<z<1 for z in fractions):
        parser.error('all fractions must lie strictly between 0 and 1')
    validate(fractions,output_path=args.output)
