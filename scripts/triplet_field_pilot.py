"""Triplet solenoid + ISR-like pair boost versus weak decay recoil, paired seeds.

Numerically integrates the published straight-string Nambu-Goto endpoint-force
approximation, or the previous central-force proxy, in a uniform Bz field.
Both endpoints are produced at the same point with zero initial internal L.
A nonzero longitudinal pair CM boost is used as a controlled ISR proxy.

This is a PHYSICS DIAGNOSTIC, not a full worldsheet or an annihilation model.
"""
import argparse
import csv
import json
import math
import statistics
import time
from dataclasses import replace
from pathlib import Path

from quirksim.field_dynamics import FieldParameters, simulate


def median(values):
    return statistics.median(values) if values else None


def run(n, params, output_dir):
    out=Path(output_dir)
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    started=time.monotonic()
    for j in range(n):
        seed=params.seed+j
        a=simulate(replace(params,seed=seed,mode='no_recoil'))
        b=simulate(replace(params,seed=seed,mode='recoil'))
        if a['plan']!=b['plan']:
            raise RuntimeError('unmatched decay plans')
        ac=[x['channel'] for x in a['decays']]
        bc=[x['channel'] for x in b['decays']]
        if ac!=bc:
            raise RuntimeError('different channels in paired histories')
        ra=a['closest_post_decay_turning_point_cm']
        rb=b['closest_post_decay_turning_point_cm']
        def value(x):
            return None if x is None else x*1e4
        firstkick=(math.sqrt(sum(z*z for z in b['decays'][0]['delta_L_from_decay_hbar']))
                   if b['decays'] else None)
        pre_mag=(math.sqrt(sum(z*z for z in b['first_decay_L_over_hbar']))
                 if b['first_decay_L_over_hbar'] is not None else None)
        rows.append({'seed':seed,'n_decays':len(bc),'channels':'+'.join(bc),
                     'pre_first_decay_L_over_hbar':pre_mag,
                     'max_pre_decay_L_over_hbar':b['max_pre_first_decay_L_over_hbar'],
                     'recoil_first_kick_abs_delta_L_over_hbar':firstkick,
                     'no_recoil_rmin_um':value(ra),'recoil_rmin_um':value(rb),
                     'paired_rmin_increase_um':None if ra is None or rb is None else (rb-ra)*1e4,
                     'no_recoil_geometric_contact':int(a['post_decay_core_encounters']>0),
                     'recoil_geometric_contact':int(b['post_decay_core_encounters']>0)})
        if (j+1)%10==0 or j+1==n:
            print(f'[field-pilot] {j+1}/{n} b={params.b_tesla}T '
                  f'force={params.force_model} elapsed={time.monotonic()-started:.1f}s',
                  flush=True)
    with (out/'paired.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    defined=[r for r in rows if r['paired_rmin_increase_um'] is not None]
    pre=[r['pre_first_decay_L_over_hbar'] for r in rows
         if r['pre_first_decay_L_over_hbar'] is not None]
    kicks=[r['recoil_first_kick_abs_delta_L_over_hbar'] for r in rows
           if r['recoil_first_kick_abs_delta_L_over_hbar'] is not None]
    result={
        'parameters':vars(params),'n':n,
        'number_with_both_closest_approaches':len(defined),
        'median_pre_decay_relative_L_over_hbar':median(pre),
        'median_first_recoil_kick_L_over_hbar':median(kicks),
        'median_no_recoil_rmin_um':median([r['no_recoil_rmin_um'] for r in defined]),
        'median_recoil_rmin_um':median([r['recoil_rmin_um'] for r in defined]),
        'median_paired_rmin_increase_um':
            median([r['paired_rmin_increase_um'] for r in defined]),
        'fraction_recoil_rmin_larger_conditional':
            sum(r['paired_rmin_increase_um']>0 for r in defined)/len(defined) if defined else None,
        'fraction_no_recoil_geometric_contacts':
            sum(r['no_recoil_geometric_contact'] for r in rows)/n,
        'fraction_recoil_geometric_contacts':
            sum(r['recoil_geometric_contact'] for r in rows)/n,
        'scope':'uniform B, same-point production, fixed ISR-like boost; '
                'straight-string approximate endpoint force, no full worldsheet; '
                'geometric core contact is not physical annihilation'}
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('[field-result] '+json.dumps(result),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',choices=['doublet','triplet'],default='triplet')
    parser.add_argument('--mass',type=float,default=225)
    parser.add_argument('--lambda-ev',type=float,default=774)
    parser.add_argument('--b-tesla',type=float,default=2)
    parser.add_argument('--beta-rel',type=float,default=.4)
    parser.add_argument('--beta-cm-x',type=float,default=.25)
    parser.add_argument('--force-model',choices=['central_proxy','ng_endpoint'],default='ng_endpoint')
    parser.add_argument('--nevents',type=int,default=50)
    parser.add_argument('--seed',type=int,default=10)
    parser.add_argument('--duration-ps',type=float,default=1000)
    parser.add_argument('--output-dir',default='results/triplet_field_pilot')
    a=parser.parse_args()
    if a.nevents<=0:
        parser.error('--nevents must be positive')
    params=FieldParameters(model=a.model,mass_gev=a.mass,lambda_ev=a.lambda_ev,
                           beta_rel=a.beta_rel,beta_cm_x=a.beta_cm_x,
                           b_tesla=a.b_tesla,force_model=a.force_model,
                           duration_ps=a.duration_ps,seed=a.seed)
    run(a.nevents,params,a.output_dir)


if __name__=='__main__':
    main()
