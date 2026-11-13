"""Mass-shell roundoff counterexample and bounded AF rerun contract."""
import json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from af_physical_weak_ensemble import prepare
from audit_physical_weak import summarize


def test_repeated_projection_and_geant_style_roundtrip(tmp_path):
    exe=tmp_path/'shell'
    subprocess.run(['g++','-O2','-std=c++17','-Wall','-Wextra','-Werror','-pedantic',str(ROOT/'research/quirkpheno_recoil/native_shell_check.cc'),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
    for value,expected in (('0','0'),('1','1'),('true','0')):
        env=dict(os.environ,QP_SPATIAL_MASS_SHELL=value)
        assert subprocess.check_output([str(exe),'enabled'],env=env,text=True).strip()==expected


def test_preflight_only_does_not_schedule_population(tmp_path):
    m=prepare(tmp_path/'probe',preflight_only=True)
    assert len(m['cases'])==2 and m['sampled_seeds']==[1] and m['spatial_mass_shell']
    assert {c['level'] for c in m['cases']}=={7,8}
    assert all(c['mode']=='recoil' and c['field_tesla']==0 and c['nominal_stop_ps']==1000 for c in m['cases'])
    assert len(list((tmp_path/'probe').glob('*.mac')))==2
    rows=[dict(name=c['name'],seed=1,field_tesla=0.,level=c['level'],mode='recoil',status='UNKNOWN_INCOMPLETE') for c in m['cases']]
    r=summarize(m,rows)
    assert r['all_transport_and_conservation_pass'] is False
    assert all(p['status']=='NOT_RUN_REQUIRED_CONTROLS' for p in r['paired_comparisons'])


def test_audit_still_rejects_boost_drift(tmp_path):
    from audit_native_noether import audit_stream
    from test_native_noether import fixture,write
    rows=fixture();rows[-1]['mechanical_j_mev_mm']=[0.]*16
    rows[-1]['mechanical_j_mev_mm'][3]=.03;rows[-1]['mechanical_j_mev_mm'][12]=-.03
    r=audit_stream(write(tmp_path,rows),0,False,150)
    assert r['status']=='FAIL' and 'orbital_mev_mm' in r['failed']
    assert r['orbital_component_max_mev_mm']['M03']==.03
    assert r['orbital_component_max_mev_mm']['M12']==0.
