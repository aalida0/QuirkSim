"""Reproducible paired stable/no-recoil/recoil geometry comparisons.

Shared seeds reuse endpoint proper-time thresholds, channel choices and decay
rest-frame randomness. Core contact is geometric, NEVER annihilation.
"""
import argparse
import csv
import json
import statistics
import time
from dataclasses import replace
from pathlib import Path
from quirksim.dynamics import Parameters, propagate


def run(n, base, output, progress_every=10, checkpoint_every=25, skip_stable=False):
    output=Path(output)
    output.mkdir(parents=True,exist_ok=True)
    rows=[]
    pair_rows=[]
    modes=('no_recoil', 'recoil') if skip_stable else ('stable', 'no_recoil', 'recoil')
    started=time.monotonic()

    def save_snapshot():
        # Atomic replacement prevents leaving truncated CSV files on interruption.
        if not rows:
            return
        for filename, values in [('events.csv',rows), ('pairs.csv',pair_rows)]:
            tmp=output/(filename+'.tmp')
            with tmp.open('w', newline='') as f:
                w=csv.DictWriter(f, fieldnames=list(values[0]))
                w.writeheader()
                w.writerows(values)
            tmp.replace(output/filename)

    for index in range(n):
        seed=base.seed+index
        event_results={}
        event_rows=[]
        for mode in modes:
            a=propagate(replace(base,seed=seed,mode=mode,store_trace=False))
            event_results[mode]=a
            channels=[e['channel'] for e in a['events']]
            rmin=a['closest_post_decay_turning_point_cm']
            event_rows.append(dict(model=base.representation,mass_gev=base.mass_neutral_gev,
                lambda_ev=base.lambda_ev,seed=seed,mode=mode,n_decays=len(channels),
                n_pion=channels.count('pion'),n_e=channels.count('electron_neutrino'),
                n_mu=channels.count('muon_neutrino'),first_decay_ps=a['first_decay_time_ps'],
                n_post_decay_turning_points=a['post_decay_encounter_count'],
                n_post_decay_core_encounters=a['post_decay_geometric_core_encounters'],
                closest_post_decay_um=None if rmin is None else 1e4*rmin,
                n_all_turning_points=len(a['turning_points']),
                n_all_core_encounters=sum(x['separation_cm']<=a['core_radius_cm']
                                         for x in a['turning_points']),
                max_separation_cm=a['maximum_separation_cm']))
        no=event_results['no_recoil']; recoil=event_results['recoil']
        if no['plan']!=recoil['plan']:
            raise AssertionError('paired event plans differ')
        pair_rows.append({'seed':seed,'n_decays':len(recoil['events']),
            'channels_equal':str([e['channel'] for e in no['events']]==
                                   [e['channel'] for e in recoil['events']]),
            'nonrecoil_closest_um':None if no['closest_post_decay_turning_point_cm'] is None
                else no['closest_post_decay_turning_point_cm']*1e4,
            'recoil_closest_um':None if recoil['closest_post_decay_turning_point_cm'] is None
                else recoil['closest_post_decay_turning_point_cm']*1e4,
            'nonrecoil_core_encounters':no['post_decay_geometric_core_encounters'],
            'recoil_core_encounters':recoil['post_decay_geometric_core_encounters']})
        rows.extend(event_rows)
        completed=index+1
        if checkpoint_every>0 and completed % checkpoint_every==0:
            save_snapshot()
        if (progress_every>0 and completed % progress_every==0) or completed==n:
            elapsed=time.monotonic()-started
            sec_per_event=elapsed/completed
            print(f'[QuirkSim] {completed}/{n} matched events ({100*completed/n:.1f}%), '
                  f'elapsed={elapsed:.1f}s, avg={sec_per_event:.2f}s/event, '
                  f'ETA≈{(n-completed)*sec_per_event/60:.1f}min',
                  flush=True)
    save_snapshot()
    summary={'benchmark': {'model':base.representation,'mass_gev':base.mass_neutral_gev,
                           'lambda_ev':base.lambda_ev,'beta':base.beta,
                           'duration_ps':base.duration_ps,'n':n},
             'mode_summaries':{},'paired':{},'scientific_scope':
             'central string proxy, conditional on initially charged pair; post-decay geometric core contact is NOT annihilation; finite observation window'}
    for mode in modes:
        rr=[r for r in rows if r['mode']==mode]
        dists=[r['closest_post_decay_um'] for r in rr if r['closest_post_decay_um'] is not None]
        summary['mode_summaries'][mode]={
            'fraction_at_least_one_decay':sum(r['n_decays']>=1 for r in rr)/n,
            'fraction_both_endpoints_decayed':sum(r['n_decays']==2 for r in rr)/n,
            'fraction_at_least_one_post_decay_geometric_core_encounter':
                sum(r['n_post_decay_core_encounters']>=1 for r in rr)/n,
            'fraction_any_geometric_core_encounter_through_end':
                sum(r['n_all_core_encounters']>=1 for r in rr)/n,
            'post_decay_closest_approach_defined_fraction':len(dists)/n,
            'median_post_decay_closest_um':statistics.median(dists) if dists else None}
    pairs=[r for r in pair_rows if r['nonrecoil_closest_um'] is not None and r['recoil_closest_um'] is not None]
    summary['paired']={
       'matching_plan_fraction':1.0,
       'fraction_matching_channels':sum(r['channels_equal']=='True' for r in pair_rows)/n,
       'both_closest_approaches_defined':len(pairs),
       'fraction_norecoil_core_but_recoil_not':sum(r['nonrecoil_core_encounters']>0 and
             r['recoil_core_encounters']==0 for r in pair_rows)/n,
       'median_recoil_minus_norecoil_um':(statistics.median(
           r['recoil_closest_um']-r['nonrecoil_closest_um'] for r in pairs)
           if pairs else None)}
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',choices=['doublet','triplet'],required=True)
    p.add_argument('--mass',type=float,required=True)
    p.add_argument('--lambda-ev',type=float,required=True)
    p.add_argument('--nevents',type=int,default=100)
    p.add_argument('--seed',type=int,default=1)
    p.add_argument('--duration-ps',type=float,default=1000)
    p.add_argument('--beta',type=float,default=.4)
    p.add_argument('--impact-cm',type=float,default=0)
    p.add_argument('--output-dir',default='results/paired')
    p.add_argument('--progress-every',type=int,default=10,
                   help='print progress and estimated remaining time every N seeds')
    p.add_argument('--checkpoint-every',type=int,default=25,
                   help='atomically write partial CSV files every N complete seeds')
    p.add_argument('--skip-stable',action='store_true',
                   help='omit redundant stable reference, run only matched no_recoil/recoil')
    a=p.parse_args()
    if a.nevents<=0: p.error('--nevents must be positive')
    base=Parameters(representation=a.model,mass_neutral_gev=a.mass,
        lambda_ev=a.lambda_ev,duration_ps=a.duration_ps,beta=a.beta,
        impact_cm=a.impact_cm,seed=a.seed)
    run(a.nevents,base,a.output_dir,a.progress_every,a.checkpoint_every,a.skip_stable)


if __name__=='__main__':
    main()
