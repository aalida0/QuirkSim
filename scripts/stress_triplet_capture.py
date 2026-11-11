"""Declared prescription sensitivity; broad gates are stress tests, not rates."""
import argparse
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.repeated_encounters import EncounterPrescription,evolve
from compare_repeated_encounters import probabilities


def calculate(source):
    original=json.loads(gzip.decompress(source.read_bytes()))
    fits=[f for f in original['fits'] if f['level']==8]
    base=original['runs'][0]['path']['model'];rows=[];budgets=[]
    for f in fits:
        e=f['fitted_energy_gev'];ell=f['fitted_L_over_hbar'];lam=base['lambda_ev']*1e-9
        budgets.append({'case':f['case'],'gate_multiple_to_open_initial_orbit':
            f['central_inner_turning_point_um']/(1.973269804e-4/base['lambda_ev']*1000),
            'minimum_directed_L_loss_per_emission_to_reach_capture':
            max(0.,ell-base['capture_max_L'])*lam/(e-base['capture_energy_in_lambda']*lam)})
        cards=[]
        for gate in (1.,10.,1e3,1e5):
            for loss in sorted({0.,1.,gate}):
                for eps in (0.,.01,.1,1.):
                    cards.append(('radiation',replace(EncounterPrescription(**base),gate_radius_multiple=gate,
                        angular_loss_per_effective_emission=loss,epsilon_per_full_period=eps)))
        for gate in (1.,1e5):
            for loss in sorted({1.,gate}):
                for cap in (1.,32.):
                    for capL in (0.,10.):
                        cards.append(('capture',replace(EncounterPrescription(**base),gate_radius_multiple=gate,
                            angular_loss_per_effective_emission=loss,epsilon_per_full_period=.1,
                            capture_energy_in_lambda=cap,capture_max_L=capL)))
        for family,model in cards:
            path=evolve(e,ell,1e13,model)
            closures=[dict(probabilities(path,spin,h,model.singlet_annihilation_width_gev),
                           singlet_fraction_at_capture=spin,excited_optical_depth_per_eligible_encounter=h)
                      for spin in (0.,.25,1.) for h in (0.,.01)]
            rows.append({'case':f['case'],'point':f['point'],'mode':f['mode'],'seed':f['seed'],
                'family':family,'gate':model.gate_radius_multiple,'loss':model.angular_loss_per_effective_emission,
                'epsilon':model.epsilon_per_full_period,'capture_energy_lambda':model.capture_energy_in_lambda,
                'capture_max_L':model.capture_max_L,'status':path['status'],
                'capture_time_ps':path['capture_time_ps'],'valid_through_ps':path['probability_valid_through_ps'],
                'energy_residual_gev':path['energy_ledger_residual_gev'],
                'angular_residual_hbar':path['angular_ledger_residual_hbar'],'closures':closures})
    paired=[];lookup={}
    for r in rows:
        key=(r['point'],r['seed'],r['family'],r['gate'],r['loss'],r['epsilon'],r['capture_energy_lambda'],r['capture_max_L'])
        lookup.setdefault(key,{})[r['mode']]=r
    for key,modes in lookup.items():
        for mode in ('collinear_forward','collinear_backward'):
            if 'isotropic' not in modes or mode not in modes:continue
            a=modes['isotropic'];b=modes[mode]
            for ac,bc in zip(a['closures'],b['closures']):
                valid=ac['status']==bc['status']=='DECLARED_MODEL_ONLY'
                paired.append({'point':a['point'],'seed':a['seed'],'control':mode,'family':a['family'],
                    'gate':a['gate'],'loss':a['loss'],'epsilon':a['epsilon'],
                    'capture_energy_lambda':a['capture_energy_lambda'],'capture_max_L':a['capture_max_L'],
                    'spin':ac.get('singlet_fraction_at_capture'),'hazard':ac.get('excited_optical_depth_per_eligible_encounter'),
                    'valid':valid,'isotropic_minus_control':ac['annihilation_probability']-bc['annihilation_probability'] if valid else None})
    return {'scope':'PRESCRIPTION_STRESS_TEST_NOT_PHYSICAL_UNCERTAINTY_BAND',
        'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'archive_sha256':original['archive_sha256'],
        'horizon_ps':1e13,'budgets':budgets,'runs':rows,'paired_contrasts':paired,
        'limitations':['Large gates and maximal directed losses are diagnostic interventions without matched emission rates.',
            'Excited hazard is per eligible encounter and remains gated.',
            'No off-gate radiation kernel, quantum spin evolution or NG backreaction is established.',
            'Truncated central paths have unknown outcome, not zero annihilation.',
            'Only the both-neutral level-8 subset is analyzed; scenarios are not statistical errors.']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=calculate(a.source);data=(json.dumps(r,separators=(',',':'),allow_nan=False)+'\n').encode()
    a.output.write_bytes(gzip.compress(data,mtime=0));print('paths',len(r['runs']),'paired scenarios',len(r['paired_contrasts']))
