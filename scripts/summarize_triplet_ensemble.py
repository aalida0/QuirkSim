"""Compact reviewer-facing summary of the audited AF kinematic-grid pilot."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import numpy as np


def summarize(archive,bounds_file,models_file):
    bounds=json.loads(bounds_file.read_text());models=json.loads(models_file.read_text())
    digest=hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest!=models['archive_sha256'] or digest!=bounds['archive_sha256']:
        raise ValueError('analysis inputs refer to different archives')
    with tarfile.open(archive) as tar:
        manifest=json.load(tar.extractfile('manifest.json'));local=json.load(tar.extractfile('local_audit.json'))
        if not local['local_checks_pass']:raise ValueError('AF local audit did not pass')
        max_four=0.;max_shell=0.;max_clock=0.;max_direction=0.;decays=0;handoffs=0;controlled=0
        for case in manifest['cases']:
            if hashlib.sha256(tar.extractfile(case['name']+'.mac').read()).hexdigest()!=case['macro_sha256']:
                raise ValueError('macro hash mismatch')
            if case['weak_enabled']:
                plan=Path(case['plan_file']).name
                if hashlib.sha256(tar.extractfile(plan).read()).hexdigest()!=case['plan_sha256']:raise ValueError('plan hash mismatch')
            records=[json.loads(line.split('QP_WEAK_JSON ',1)[1]) for line in tar.extractfile(case['name']+'.log').read().decode().splitlines() if 'QP_WEAK_JSON ' in line]
            for r in records:
                if r['kind']=='handoff':
                    if not r['history_preserved']:raise ValueError('lost handoff history')
                    handoffs+=1;continue
                decays+=1;p=np.array(r['parent_mev']);d=np.array(r['daughter_mev']);pi=np.array(r['pion_mev'])
                max_four=max(max_four,float(np.max(abs(p-d-pi))))
                max_clock=max(max_clock,abs(r['proper_ps']-r['threshold_ps']))
                for label,v in [('parent',p),('daughter',d),('pion',pi)]:
                    max_shell=max(max_shell,abs(v[0]-np.sqrt(r[label+'_mass_mev']**2+v[1:]@v[1:])))
                if case.get('direction_mode')=='collinear_cm':
                    bz=case['pair_boost_z'];gamma=1/np.sqrt(1-bz*bz);pc=p.copy();kc=pi.copy()
                    pc[0]=gamma*(p[0]-bz*p[3]);pc[3]=gamma*(p[3]-bz*p[0])
                    kc[0]=gamma*(pi[0]-bz*pi[3]);kc[3]=gamma*(pi[3]-bz*pi[0])
                    b=pc[1:]/pc[0];b2=b@b;gg=pc[0]/r['parent_mass_mev'];rest=kc[1:].copy()
                    if b2:rest+=((gg-1)*(b@kc[1:])/b2-gg*kc[0])*b
                    norm=np.linalg.norm(pc[1:]);q=0 if r['parent_pdg']>0 else 1
                    target=(pc[1:]/norm if norm>1e-14*r['parent_mass_mev'] else np.array(case['plan'][q]['direction']))*case['collinear_sign']
                    max_direction=max(max_direction,float(np.linalg.norm(rest/np.linalg.norm(rest)-target)));controlled+=1
        if decays!=handoffs or max_four>1e-6 or max_clock>2e-9 or max_direction>1e-8:raise ValueError('independent logged-vertex checks failed')
    by_name={c['name']:c for c in bounds['cases']};points=[]
    for point,counts in manifest['ensemble_summary'].items():
        modes={}
        for mode in ('isotropic','collinear_forward','collinear_backward'):
            fits=[f for f in models['fits'] if f['point']==point and f['mode']==mode and f['level']==8]
            minima=np.array([f['separation_um'] for f in fits]);status_counts={}
            for fit in fits:
                status=by_name[fit['case']]['status'];status_counts[status]=status_counts.get(status,0)+1
            modes[mode]={'neutral_qualified_seeds':[f['seed'] for f in fits],
                         'CM_chord_minimum_um':{'min':float(minima.min()),'median':float(np.median(minima)),'max':float(minima.max())} if len(minima) else None,
                         'lab_speed_bound_classifications':status_counts}
        selected=[s for s in models['paired_conditional_direction_comparisons'] if s['point']==point and
            s['epsilon']==.1 and s['angular_loss']==1. and s['horizon_ps']==1e13 and
            s['capture_singlet_fraction']==.25 and s['excited_optical_depth']==0.]
        contrasts=[s['mean_isotropic_minus_control_probability'] for s in models['paired_conditional_direction_comparisons'] if s['point']==point]
        points.append({'point':point,'isotropic_population_counts':counts,'modes':modes,
                       'selected_declared_model_contrasts':selected,
                       'mean_paired_contrast_range_over_tested_closures':[min(contrasts),max(contrasts)] if contrasts else None})
    lookup={(f['point'],f['mode'],f['seed'],f['level']):f for f in models['fits']}
    shifts=[]
    for (point,mode,seed,level),fit in lookup.items():
        if level!=7:continue
        fine=lookup.get((point,mode,seed,8))
        if fine:shifts.append({'point':point,'mode':mode,'seed':seed,
                              'CM_minimum_shift_nm':1000*(fine['separation_um']-fit['separation_um'])})
    projections=[f['angular_fit_diagnostic'] for f in models['fits'] if f['angular_fit_diagnostic']['symmetry_projection_applied']]
    return {'scope':'KINEMATIC_GRID_PILOT_CONDITIONAL_NEUTRAL_SUBSET',
        'archive_sha256':digest,'transport_cases':len(manifest['cases']),
        'AF_build_passed':manifest['build']['returncode']==0,
        'AF_local_audit_passed':local['local_checks_pass'],
        'independently_checked_logged_decays':decays,'preserved_handoffs':handoffs,
        'max_four_component_residual_mev':max_four,'max_shell_energy_residual_mev':max_shell,
        'max_proper_threshold_residual_ps':max_clock,'points':points,
        'independently_checked_control_decays':controlled,'max_control_direction_error':max_direction,
        'refinement_checks':shifts,
        'angular_floor_correction':{'projected_control_windows':len(projections),
            'max_transverse_position_residual_mm':max(p['transverse_position_residual_mm'] for p in projections) if projections else None,
            'scope':'exact zero-field CM-collinear symmetry only; raw residuals retained in model artifact'},
        'declared_orbit_scenarios':len(models['runs']),
        'interpretation':'2T controls already miss the gate in the central model; zero-field directional comparison shows the transverse-kick mechanism',
        'physical_inclusive_annihilation_probability':None,
        'limitations':['8 shared seeds across four unweighted points; only 3 or 5 neutral-qualified seeds per point',
            'forced pion channel and conditional angular controls',
            '2T control nanometre-scale exclusion is a convergent interpolation/model result; speed-only bounds remain unresolved',
            'capture and spin closure choices can remove the annihilation contrast',
            'neutral subset is not an inclusive population survival estimate']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path)
    p.add_argument('--bounds',type=Path,required=True);p.add_argument('--models',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=summarize(a.archive,a.bounds,a.models);a.output.write_text(json.dumps(r,indent=2)+'\n')
    print('AF cases',r['transport_cases'],'checked decays',r['independently_checked_logged_decays'],'orbit scenarios',r['declared_orbit_scenarios'])
