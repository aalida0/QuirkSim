"""Disentangle integrator/short-distance-regulator error from NG mesh error.

All comparisons condition on the SAME total spatial arclength crossing,
not on fixed per-element distances. This is a pre-pinch finite-element
study; no annihilation, worldsheet reconnection, or radiation is modeled.
This script reports a PHYSICS BLOCK when endpoint-chord uncertainty
is large relative to hbar*c/Lambda_IC. CI failure is reserved for numerical
inconsistency, not for an honestly reported physics limitation.
"""
import argparse
import json
from pathlib import Path
import math
from equidistributed_gauge_stress import case
from quirksim.dynamics import HBARC_GEV_CM

REFERENCE_CORE_CM=HBARC_GEV_CM/(774e-9)


def run(cutoff_fraction=2e-5,nodes=(5,7,11,15),
        horizon_ps=1000,output='results/near_core_error_budget.json'):
    numerical=[]
    for n in nodes:
        settings=[
            ('baseline',0.5,2e-9,1e-10,1e-9),
            ('tight_step',0.125,2e-10,1e-11,1e-9),
            ('half_regulator',0.5,2e-9,1e-10,5e-10)
        ]
        row={'nodes':n,'settings':{}}
        for name,step,rtol,atol,reg in settings:
            result=case(n,2.,.25,horizon_ps,pinch_fraction=cutoff_fraction,
                        max_step_ps=step,rtol=rtol,atol=atol,regulator_cm=reg)
            row['settings'][name]=result
            print('[near-core-case] '+json.dumps({'nodes':n,'setting':name,
                'status':result['status'],
                'time_ps':result.get('time_ps'),
                'chord_cm':result.get('endpoint_separation_cm'),
                'total_length_cm':result.get('total_string_polyline_length_cm'),
                'energy_drift_gev':result.get('energy_drift_gev'),
                'reason':result.get('reason')}),flush=True)
        for name,r in row['settings'].items():
            if r['status']!='PHYSICAL_ARCLENGTH_CUTOFF' or r['stops']!=['total_arclength_cutoff']:
                raise AssertionError('numerical integration not valid: n='+str(n)+' '+name+
                                     ' '+r['status'])
        b=row['settings']['baseline']
        tight=row['settings']['tight_step']
        reg=row['settings']['half_regulator']
        row['absolute_chord_step_error_cm']=abs(
            tight['endpoint_separation_cm']-b['endpoint_separation_cm'])
        row['absolute_chord_regulator_shift_cm']=abs(
            reg['endpoint_separation_cm']-b['endpoint_separation_cm'])
        row['absolute_stop_time_step_error_ps']=abs(tight['time_ps']-b['time_ps'])
        row['absolute_stop_time_regulator_shift_ps']=abs(reg['time_ps']-b['time_ps'])
        row['max_energy_drift_gev']=max(abs(r['energy_drift_gev']) for r in row['settings'].values())
        numerical.append(row)
        print('[near-core-errors] '+json.dumps({
            k:v for k,v in row.items() if k!='settings'}),flush=True)

    chords={n['nodes']:n['settings']['tight_step']['endpoint_separation_cm']
            for n in numerical}
    mesh_spread=max(chords.values())-min(chords.values())
    errors={
       'max_timestep_chord_shift_cm':max(r['absolute_chord_step_error_cm'] for r in numerical),
       'max_regulator_chord_shift_cm':max(r['absolute_chord_regulator_shift_cm'] for r in numerical),
       'mesh_chord_spread_cm':mesh_spread,
       'reference_1_over_lambda_cm':REFERENCE_CORE_CM,
       'mesh_spread_in_reference_units':mesh_spread/REFERENCE_CORE_CM,
       'max_timestep_shift_in_reference_units':max(r['absolute_chord_step_error_cm']
                                                    for r in numerical)/REFERENCE_CORE_CM,
       'max_regulator_shift_in_reference_units':max(r['absolute_chord_regulator_shift_cm']
                                                     for r in numerical)/REFERENCE_CORE_CM
    }
    # Endpoints must converge independently of the arbitrary numerical
    # integration settings; preserving total energy alone is insufficient.
    numerical_pass=(errors['max_timestep_chord_shift_cm']<0.1*REFERENCE_CORE_CM
                    and errors['max_regulator_chord_shift_cm']<0.1*REFERENCE_CORE_CM)
    core_pass=(numerical_pass and len(nodes)>=3 and
               mesh_spread<0.1*REFERENCE_CORE_CM)
    report={
       'status':('CORE_REFERENCE_GEOMETRY_CONVERGED_DIAGNOSTIC_ONLY' if core_pass
                  else 'SPATIAL_CONVERGENCE_NOT_ASSESSED' if len(nodes)<3
                  else 'CORE_REFERENCE_GEOMETRY_BLOCKED'),
       'numerical_stability_gate':bool(numerical_pass),
       'core_reference_geometry_gate':bool(core_pass),
       'cutoff_fraction_of_initial_length':cutoff_fraction,
       'total_arclength_cutoff_cm':.1*cutoff_fraction,
       'chord_by_nodes_tight_step_cm':chords,
       'error_budget':errors,'samples':numerical,
       'limitations':[
          'core reference hbar*c/Lambda_IC is not physical annihilation radius',
          'pre-stretched string, fixed 2T and boost slice, no weak decays',
          'no dynamical string continuation past first cutoff',
          'midpoint Nambu-Goto finite elements, not continuum worldsheet certified',
          'no contact/annihilation probability predicted']
    }
    dest=Path(output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,indent=2)+'\n')
    print('[near-core-summary] '+json.dumps({
       'status':report['status'],
       'numerical_stability_gate':numerical_pass,
       'core_reference_geometry_gate':core_pass,
       'chords':chords,'error_budget':errors}),flush=True)
    if not numerical_pass:
        raise AssertionError('numerical timestep/regulator core reference gate failed')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--nodes',default='5,7,11,15')
    p.add_argument('--cutoff-fraction',type=float,default=2e-5)
    p.add_argument('--output',default='results/near_core_error_budget.json')
    args=p.parse_args()
    ns=tuple(int(v.strip()) for v in args.nodes.split(','))
    if not ns or min(ns)<3 or max(ns)>32 or len(set(ns))!=len(ns):
        p.error('nodes must be distinct between 3 and 32')
    if not 0<args.cutoff_fraction<1:
        p.error('invalid cutoff fraction')
    run(args.cutoff_fraction,ns,output=args.output)
