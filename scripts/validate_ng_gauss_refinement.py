"""Validate independently produced Gauss2/Gauss3 NG refinement results.

Inputs: run scripts/equidistributed_gauge_stress.py at common parameters
for N=13,17,21, quadrature=2, and N=21, quadrature=3.

The <=0.1*(hbar*c/Lambda_IC) criterion is a *numerical diagnostic* only.
Passing it is NOT a validation of the continuum string through contact,
the quirk annihilation rate, or any collider acceptance.
"""
import argparse
import json
import math
from pathlib import Path
from quirksim.dynamics import HBARC_GEV_CM

CORE_REF_CM=HBARC_GEV_CM/(774e-9)
EXPECTED_ARCLENGTH_CM=.1*2e-5


def validate(input_dir,output=None):
    root=Path(input_dir)
    cases={}
    for n,q in ((13,2),(17,2),(21,2),(21,3)):
        path=root/f'ng_n{n}_q{q}.json'
        data=json.loads(path.read_text())
        if len(data.get('records',[]))!=1:
            raise AssertionError(f'{path}: expected exactly one benchmark')
        row=data['records'][0]
        if (row.get('nodes')!=n or row.get('spatial_quadrature_order')!=q
                or row.get('status')!='PHYSICAL_ARCLENGTH_CUTOFF'
                or row.get('stops')!=['total_arclength_cutoff']):
            raise AssertionError(f'{path}: unexpected solver status/identity')
        if not (math.isfinite(row['time_ps']) and row['time_ps']<1000
                and row['b_tesla']==2.0 and row['beta_boost_x']==.25
                and row['gauge_point_mass_gev']==0):
            raise AssertionError(f'{path}: wrong benchmark or gauge')
        if (abs(row['numerical_regulator_cm']-2e-10)>1e-20
                or abs(row['common_total_length_pinch_cutoff_cm']-
                       EXPECTED_ARCLENGTH_CM)>1e-15
                or abs(row['total_string_polyline_length_cm']-
                       EXPECTED_ARCLENGTH_CM)>2e-12):
            raise AssertionError(f'{path}: mismatched cutoff or regulator')
        if abs(row['energy_drift_gev'])>1e-7:
            raise AssertionError(f'{path}: excessive energy drift')
        if row['max_equal_length_constraint_cm']>1e-6:
            raise AssertionError(f'{path}: equal-arclength gauge not maintained')
        cases[(n,q)]=row

    spatial={n:cases[(n,2)]['endpoint_separation_cm'] for n in (13,17,21)}
    chord_spread=max(spatial.values())-min(spatial.values())
    times=[cases[(n,2)]['time_ps'] for n in (13,17,21)]
    quad_chord_error=abs(cases[(21,2)]['endpoint_separation_cm']-
                         cases[(21,3)]['endpoint_separation_cm'])
    quad_time_error=abs(cases[(21,2)]['time_ps']-
                        cases[(21,3)]['time_ps'])
    spatial_pass=chord_spread < .1*CORE_REF_CM
    quad_pass=quad_chord_error < .01*CORE_REF_CM
    time_pass=max(times)-min(times)<5e-4 and quad_time_error<5e-4
    report={
        'status':('SELECTED_PRE_CONTACT_GEOMETRY_GATE_PASS'
                  if spatial_pass and quad_pass and time_pass else
                  'BLOCKED_SELECTED_PRE_CONTACT_GEOMETRY'),
        'mesh_vertex_counts':[13,17,21],
        'spatial_quadrature_orders':[2,3],
        'selected_triplet_mass_gev':225,
        'selected_lambda_ev':774,
        'field_tesla':2,
        'pair_boost_beta_x':.25,
        'stop_arclength_cm':EXPECTED_ARCLENGTH_CM,
        'regulator_cm':2e-10,
        'reference_hbarc_over_lambda_cm':CORE_REF_CM,
        'chord_by_nodes_cm':spatial,
        'max_mesh_chord_spread_cm':chord_spread,
        'mesh_chord_spread_over_reference':chord_spread/CORE_REF_CM,
        'gauss2_vs_gauss3_21node_chord_difference_cm':quad_chord_error,
        'gauss2_vs_gauss3_chord_difference_over_reference':quad_chord_error/CORE_REF_CM,
        'max_mesh_hitting_time_spread_ps':max(times)-min(times),
        'quadrature_hitting_time_difference_ps':quad_time_error,
        'mesh_geometry_gate':spatial_pass,
        'quadrature_gate':quad_pass,
        'time_gate':time_pass,
        'worldsheet_contact_continuation_gate':'BLOCKED_NOT_IMPLEMENTED',
        'annihilation_and_radiation_gate':'BLOCKED_NOT_IMPLEMENTED',
        'scope':'selected 2T boosted pre-stretched triplet, no weak decay; '
                'pre-contact finite-element mesh and quadrature validation only'
    }
    if output:
        out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
        out.write_text(json.dumps(report,indent=2)+'\n')
    print('[ng-selected-geometry-gate] '+json.dumps(report),flush=True)
    if report['status']!='SELECTED_PRE_CONTACT_GEOMETRY_GATE_PASS':
        raise AssertionError('NG selected-slice core-reference geometric convergence BLOCKED')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--input-dir',default='results/ng-gauss-validated')
    p.add_argument('--output',default='results/ng_selected_geometry_gate.json')
    a=p.parse_args()
    validate(a.input_dir,a.output)
