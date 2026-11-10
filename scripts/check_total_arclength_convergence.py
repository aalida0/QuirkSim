"""Fail-closed 2T triplet convergence at a COMMON TOTAL string arclength.

This is an independent CLI validation gate over N=5/7/11; it must NOT
mistake successful GitHub job execution for physical 1000ps completion.
All cases evolve until the earliest near-zero-length boundary, not through
an unvalidated folded/reconnected/annihilating worldsheet.
"""
import json
import math
from pathlib import Path
from equidistributed_gauge_stress import case


def validate(cutoff_fractions=(1e-3,4e-4), horizon_ps=1000):
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
    Path('results').mkdir(exist_ok=True)
    Path('results/constrained_ng_common_arclength.json').write_text(
        json.dumps(report,indent=2)+'\n')
    print('[common-arclength-summary] '+json.dumps({
       'status':report['status'],
       'trials':[{'cutoff':t['common_total_cutoff_cm'],
                  'spread_ps':t['max_stop_time_spread_ps'],
                  'passed':t['convergence_gate_pass']} for t in trials]}),flush=True)
    if report['status']!='PRE_PINCH_CONVERGENCE_VALIDATED_ONLY':
        raise AssertionError('common-arclength mesh convergence validation failed')
    return report


if __name__=='__main__':
    validate()
