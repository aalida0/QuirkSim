"""Nonradial triplet sensitivity scan for the *central-string proxy*.

A nonzero transverse initial offset b creates initial angular momentum
L_y = -(b * p_z) / (hbar*c), with back-to-back longitudinal momenta.
This is a controlled initial-condition stress test, **not** a realistic
proton-collision production-vertex distribution.

For equal-mass particles created at the same point, giving the pair an
opening angle does not by itself generate internal orbital angular momentum
in their center-of-momentum frame.
"""
import argparse
import csv
import json
import math
import statistics
import time
from dataclasses import replace
from pathlib import Path

from quirksim.physics import benchmark
from quirksim.dynamics import Parameters, HBARC_GEV_CM, propagate


def _median(seq):
    return statistics.median(seq) if seq else None


def impact_scan(impact_um, n, base, output_dir, progress_every=10):
    if n <= 0:
        raise ValueError('nevents must be positive')
    if any(b < 0 or not math.isfinite(b) for b in impact_um):
        raise ValueError('impact offsets must be finite and nonnegative')
    out=Path(output_dir)
    out.mkdir(parents=True,exist_ok=True)
    m=benchmark(base.representation,base.mass_neutral_gev)['charged_mass_gev']
    pz=m*base.beta/math.sqrt(1-base.beta**2)
    start=time.monotonic()
    rows=[]
    events=[]
    for impact in impact_um:
        impact_cm=impact*1e-4
        deviations=[]
        norecoil=[]
        recoil=[]
        positive=0
        positive_n=0
        channel_matched=0
        recoil_core=0
        control_core=0
        both_defined=0
        any_decay=0
        for j in range(n):
            seed=base.seed+j
            params=replace(base,impact_cm=impact_cm,seed=seed)
            no=propagate(replace(params,mode='no_recoil',store_trace=False))
            yes=propagate(replace(params,mode='recoil',store_trace=False))
            if no['plan'] != yes['plan']:
                raise AssertionError('paired random plan mismatch')
            if [e['channel'] for e in no['events']] == [e['channel'] for e in yes['events']]:
                channel_matched+=1
            a=no['closest_post_decay_turning_point_cm']
            b=yes['closest_post_decay_turning_point_cm']
            any_decay+=bool(yes['events'])
            control_core+=int(no['post_decay_geometric_core_encounters']>0)
            recoil_core+=int(yes['post_decay_geometric_core_encounters']>0)
            if a is not None and b is not None:
                both_defined+=1
                d=(b-a)*1e4
                deviations.append(d)
                norecoil.append(a*1e4)
                recoil.append(b*1e4)
                positive+=int(d>0)
                positive_n+=int(d>0.1) # >0.1 micrometre increase
            events.append({
                'impact_um':impact,'seed':seed,'channel_matched':
                    [e['channel'] for e in no['events']]==[e['channel'] for e in yes['events']],
                'n_decays':len(yes['events']),'no_recoil_closest_um':None if a is None else a*1e4,
                'recoil_closest_um':None if b is None else b*1e4,
                'difference_um':None if a is None or b is None else (b-a)*1e4,
                'no_recoil_core':int(no['post_decay_geometric_core_encounters']>0),
                'recoil_core':int(yes['post_decay_geometric_core_encounters']>0)})
            if progress_every>0 and (j+1)%progress_every==0:
                elapsed=time.monotonic()-start
                print(f'[impact] b={impact:g} um: {j+1}/{n} seeds, '
                      f'elapsed total={elapsed:.1f}s',flush=True)
        rows.append({
            'impact_um':impact,'initial_L_over_hbar_abs':impact_cm*pz/HBARC_GEV_CM,
            'n':n,'fraction_with_decay':any_decay/n,
            'channel_matched_fraction':channel_matched/n,
            'both_post_decay_minima_defined':both_defined,
            'median_no_recoil_closest_um':_median(norecoil),
            'median_recoil_closest_um':_median(recoil),
            'median_paired_difference_um':_median(deviations),
            'fraction_recoil_rmin_larger_conditional':
                (positive/both_defined if both_defined else None),
            'fraction_recoil_increase_over_0p1um_conditional':
                (positive_n/both_defined if both_defined else None),
            'fraction_with_post_decay_geometric_core_encounter_no_recoil':
                control_core/n,
            'fraction_with_post_decay_geometric_core_encounter_recoil':
                recoil_core/n})
        last=rows[-1]
        print(f"[impact-result] b={impact:g} um | L/hbar={last['initial_L_over_hbar_abs']:.3e} "
              f"| med rmin no={last['median_no_recoil_closest_um']} um "
              f"recoil={last['median_recoil_closest_um']} um "
              f"paired_delta={last['median_paired_difference_um']} um "
              f"| core(no,recoil)=({control_core}/{n},{recoil_core}/{n})",flush=True)
        # Store a completed-offset checkpoint, not a partial seed sample.
        def write_csv(filename,records):
            tmp=out/(filename+'.tmp')
            with tmp.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
            tmp.replace(out/filename)
        write_csv('summary.csv',rows)
        write_csv('paired_events.csv',events)
    results={'model':'triplet' if base.representation=='triplet' else base.representation,
             'mass_gev':base.mass_neutral_gev,'lambda_ev':base.lambda_ev,
             'beta':base.beta,'duration_ps':base.duration_ps,
             'n_per_impact':n,'seed_first':base.seed,
             'core_radius_um':HBARC_GEV_CM/(base.lambda_ev*1e-9)*1e4,
             'interpretation':('controlled transverse-offset/initial-angular-momentum '
                               'sensitivity of a central-string potential only, '
                               'NOT realistic pp production, radiation or annihilation'),
             'results':rows}
    (out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--impact-um',default='0,0.01,0.1,1,3,10',
                   help='comma-separated initial transverse separations in micrometres')
    p.add_argument('--mass',type=float,default=225)
    p.add_argument('--lambda-ev',type=float,default=774)
    p.add_argument('--beta',type=float,default=0.4)
    p.add_argument('--duration-ps',type=float,default=1000)
    p.add_argument('--nevents',type=int,default=20)
    p.add_argument('--seed',type=int,default=101)
    p.add_argument('--progress-every',type=int,default=10)
    p.add_argument('--output-dir',default='results/triplet_nonradial')
    a=p.parse_args()
    offsets=[float(x.strip()) for x in a.impact_um.split(',')]
    if not offsets: p.error('empty offsets')
    base=Parameters(representation='triplet',mass_neutral_gev=a.mass,
                    lambda_ev=a.lambda_ev,beta=a.beta,duration_ps=a.duration_ps,seed=a.seed)
    result=impact_scan(offsets,a.nevents,base,a.output_dir,a.progress_every)
    print('[impact-done] '+json.dumps(result['results']),flush=True)


if __name__=='__main__':
    main()
