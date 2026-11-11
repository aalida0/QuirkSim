"""Audit saved two-resolution AF pilots without extracting or rerunning them."""
import argparse
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import tarfile


def vector_distance(a, b):
    return math.sqrt(sum((x-y)**2 for x, y in zip(a, b)))


def endpoint_difference(a, b, key):
    return max(vector_distance(x, y) for x, y in zip(a[key], b[key]))


def audit(archive):
    results, provenance = {}, {}
    with tarfile.open(archive) as handle:
        for entry in handle.getmembers():
            if not entry.isfile():
                continue
            path = PurePosixPath(entry.name)
            if path.name in ('evolution.json', 'provenance.json'):
                case = path.parts[-3]
                target = results if path.name == 'evolution.json' else provenance
                if case in target:
                    raise ValueError(f'Duplicate archived case: {case}')
                target[case] = json.load(handle.extractfile(entry))
    if len(results) != 8 or results.keys() != provenance.keys():
        raise ValueError('Expected eight distinct results with provenance')
    cases, comparisons = [], []
    for name, r in sorted(results.items()):
        # All submitted cards use the explicit 8-reference-length matching surface.
        reference = r['matching_radius_cm']/8
        first_failure = next((row['time_ps'] for row in r['records']
                              if row['chord_closure_error_cm'] >= .1*reference), None)
        cases.append(dict(case=name, source_commit=provenance[name]['source_commit'],
                          termination=r['termination'], final_time_ps=r['final']['time_ps'],
                          numerical_validation=r['numerical_validation'],
                          first_saved_closure_failure_ps=first_failure,
                          decays=[dict(endpoint=e['endpoint'], time_ps=e['time_ps']) for e in r['decays']],
                          reflected_discontinuities=len(r['wave_reflections']),
                          final_separation_cm=vector_distance(*r['final']['positions_cm']),
                          final_energy_coordinate_interval_ps=r['final']['energy_coordinate_interval_ps'],
                          energy_residual_gev=r['energy_with_weak_products_residual_gev']))
    for model in ('triplet', 'doublet'):
        for seed in (1, 25):
            a, b = [results[f'{model}_recoil_seed{seed}_level{level}'] for level in (0, 1)]
            reference = a['matching_radius_cm']/8
            pairs = []
            for x in a['records']:
                y = min(b['records'], key=lambda row: abs(row['time_ps']-x['time_ps']))
                mismatch = abs(x['time_ps']-y['time_ps'])
                if mismatch <= 1e-7:
                    pairs.append(dict(time_ps=x['time_ps'], time_mismatch_ps=mismatch,
                                      position_difference_cm=endpoint_difference(x, y, 'positions_cm'),
                                      momentum_difference_gev=endpoint_difference(x, y, 'momenta_gev')))
            worst = max(pairs, key=lambda row: row['position_difference_cm'])
            decay_pairs = list(zip(a['decays'], b['decays']))
            comparisons.append(dict(model=model, seed=seed,
                comparison_method='nearest saved times within 1e-7 ps; no interpolation',
                matched_record_count=len(pairs), last_matched_time_ps=pairs[-1]['time_ps'],
                max_time_mismatch_ps=max(row['time_mismatch_ps'] for row in pairs),
                worst_position_comparison=worst,
                max_position_difference_in_reference_lengths=worst['position_difference_cm']/reference,
                final_times_equal=a['final']['time_ps'] == b['final']['time_ps'],
                final_endpoint_difference_cm=(endpoint_difference(a['final'], b['final'], 'positions_cm')
                    if a['final']['time_ps'] == b['final']['time_ps'] else None),
                decay_endpoint_sequence_equal=[e['endpoint'] for e in a['decays']] == [e['endpoint'] for e in b['decays']],
                max_decay_time_difference_ps=max((abs(x['time_ps']-y['time_ps']) for x, y in decay_pairs), default=None),
                reflection_counts=[len(r['wave_reflections']) for r in (a, b)],
                last_saved_comparison=pairs[-1]))
    return dict(archive_sha256=hashlib.sha256(Path(archive).read_bytes()).hexdigest(),
                scope='SAVED_CLASSICAL_NUMERICAL_DIAGNOSTICS_ONLY',
                reference_definition='matching radius / 8 for these fixed submitted cards',
                cases=cases, comparisons=comparisons,
                physical_survival_probability=None,
                caveats=['Saved records do not resolve every encounter minimum.',
                         'Two resolutions are a consistency check, not an observed convergence order.',
                         'Stop reasons do not identify their triggering internal RK/event state.'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = audit(args.archive)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({k: output[k] for k in ('archive_sha256', 'scope')}))
