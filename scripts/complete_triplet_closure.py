"""Combine G-even compact S/P widths with saved capture paths and population bounds."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tarfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.physical_closure import compact_subspace_interval
from quirksim.quirk_g_parity import neutral_capture_channels
from scipy.stats import beta


def cp_interval(k,n,alpha=.05):
    if not 0<=k<=n or n<1:raise ValueError('valid binomial counts required')
    return [0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1)),
            1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k))]


def calculate(frozen,pwave,compact,archive):
    frozen_data=json.loads(frozen.read_text());p=json.loads(pwave.read_text());s=json.loads(compact.read_text())
    if p['compact_input_sha256']!=hashlib.sha256(compact.read_bytes()).hexdigest():raise ValueError('incompatible compact inputs')
    if p['model']!=s['model']:raise ValueError('different compact model identities')
    if frozen_data['archive_sha256']!=hashlib.sha256(archive.read_bytes()).hexdigest():raise ValueError('wrong archive')
    widths={'1S0':s['sum_modeled_two_boson_widths_gev'],**{'3P'+j:w for j,w in p['states'][-1]['sum_modeled_widths_gev'].items()}}
    if any(not c['allowed_from_source_with_G_even_radiation'] for name,c in neutral_capture_channels().items() if name in widths):raise ValueError('forbidden compact support')
    cases=[]
    for c in frozen_data['cases']:
        captured=c['capture_time_ps'] is not None
        age=1e13-c['capture_time_ps'] if captured else 0.
        intervals={str(eta):compact_subspace_interval(float(captured)*eta,age,widths.values()) for eta in (0.,.1,1.)}
        cases.append({'case':c['case'],'point':c['point'],'mode':c['mode'],'seed':c['seed'],
            'capture_time_ps':c['capture_time_ps'],'baseline_compact_occupation_envelopes_by_matching_fraction':intervals,
            'one_ns_capture_only_annihilation':c['one_ns_post_handoff']['annihilation_interval']})
    with tarfile.open(archive) as t:m=json.load(t.extractfile('manifest.json'))
    b=m['reference_pion_branching'];population=[]
    for point,pop in m['ensemble_summary'].items():
        n=pop['population_seeds'];qualified=[c for c in cases if c['point']==point and c['mode']=='isotropic']
        lower=sum(c['baseline_compact_occupation_envelopes_by_matching_fraction']['1.0']['lower'] for c in qualified)/n
        upper=(sum(c['baseline_compact_occupation_envelopes_by_matching_fraction']['1.0']['upper'] for c in qualified)+n-len(qualified))/n
        population.append({'point':point,'all_isotropic_seeds':n,'both_neutral_qualified':len(qualified),
            'neutral_fraction_exact_binomial_95pct_interval':cp_interval(len(qualified),n),
            'conditional_forced_pion_sample_annihilation_interval_all_categories':[lower,upper],
            'all_branching_channels_sample_interval':[b*b*lower,b*b*upper+1-b*b],
            'statistical_interpretation':'pilot seed sample only; interval is a missing-category envelope, not a confidence interval or production prediction'})
    return {'scope':'CLOSED_DECLARED_G_EVEN_LOWEST_COMPACT_SUBSPACE_NOT_INCLUSIVE_PHYSICAL_RATE',
        'inputs_sha256':{str(x.name):hashlib.sha256(x.read_bytes()).hexdigest() for x in (frozen,pwave,compact,archive)},
        'compact_widths_gev':widths,'cases':cases,'population_accounting':population,
        'pion_branching':b,'missing_non_pion_pair_probability':1-b*b,
        'assumptions':['Minimal two-color vectorlike Y=0 Lagrangian with exact G and G-even production/radiation.',
            'At matching, compact support is restricted to lowest 1S0 and 3P_J levels; higher states remain outside this closure.',
            'Matching fractions are scenarios, not measured probabilities; eta=1 is the ideal absorbing boundary.',
            'Widths are LO modeled two-boson widths, with no imported scalar/Majorana rates.',
            'Non-pion channels and 0/1-decay histories retain explicit unknown outcomes.',
            'CP intervals assume independent lifetime/direction seeds at each fixed point; points share seeds and are not independent.',
            'No PDF weights or inclusive detector-time annihilation probability are inferred.'],
        'release':{'specified_compact_spin_sector_calculation':'COMPLETE_AT_DECLARED_LO_TWO_VECTOR_ORDER',
            'conditional_capture_to_annihilation_implementation':'COMPLETE_FOR_DECLARED_SUBSPACE',
            'unique_physical_formation_probability':'NOT_ESTABLISHED',
            'production_weighted_inclusive_result':'REQUIRES_ADDITIONAL_INPUT_AND_TRANSPORT'}}


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]/'studies/data';a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('archive',type=Path);a.add_argument('--output',type=Path,required=True);v=a.parse_args()
    result=calculate(root/'triplet_frozen_closure_2026-11-12.json',root/'triplet_pwave_2026-11-12.json',root/'neutral_triplet_quantum_2026-11-11.json',v.archive)
    v.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('compact widths GeV',result['compact_widths_gev']);print('unmodeled non-pion fraction',result['missing_non_pion_pair_probability'])
