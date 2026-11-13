"""Invariant, failure-injection and campaign checks, independent of Geant4."""
import importlib.util,json,subprocess,sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from audit_native_noether import audit_stream


def test_native_characteristic_reflections(tmp_path):
    exe=tmp_path/'noether'
    subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic',str(ROOT/'research/quirkpheno_recoil/native_noether_check.cc'),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)


def fixture():
    d=dict(schema=1,reason='transport',step=1,vertices=0,endpoint_ct_mm=[50.,50.],contour_closure_mm=[0.]*4,mechanical_p_mev=[500000.,0.,0.,0.],emitted_p_mev=[0.]*4,initial_p_mev=[500000.,0.,0.,0.],mechanical_j_mev_mm=[0.]*16,emitted_j_mev_mm=[0.]*16,initial_j_mev_mm=[0.]*16,external_torque_z_mev_mm=0.)
    return [d,dict(d,step=2000),dict(d,step=4000,reason='end_event')]


def write(tmp_path,rows):
    p=tmp_path/'stream.jsonl';p.write_text(''.join(json.dumps(d)+'\n' for d in rows));return p


def test_audit_missing_and_failure_injection(tmp_path):
    rows=fixture();assert audit_stream(write(tmp_path,rows),0,False,150)['status'].startswith('PASS')
    rows[-1]['mechanical_p_mev']=[499000,0,0,0]
    assert 'energy_mev' in audit_stream(write(tmp_path,rows),0,False,150)['failed']
    with pytest.raises(ValueError):audit_stream(write(tmp_path,rows[:2]),0,False,150)
    rows=fixture();rows[-1]['mechanical_p_mev']=[float('nan'),0,0,0]
    with pytest.raises(ValueError):audit_stream(write(tmp_path,rows),0,False,150)
    rows=fixture();rows[-1]['initial_p_mev']=[500001,0,0,0]
    with pytest.raises(ValueError):audit_stream(write(tmp_path,rows),0,False,150)


def test_axial_external_torque_not_full_j_conservation(tmp_path):
    rows=fixture()
    rows[-1]['mechanical_j_mev_mm']=[0.]*16;rows[-1]['mechanical_j_mev_mm'][6]=12.;rows[-1]['mechanical_j_mev_mm'][9]=-12.
    rows[-1]['external_torque_z_mev_mm']=12.
    assert audit_stream(write(tmp_path,rows),2,False,150)['status'].startswith('PASS')
    with pytest.raises(ValueError):audit_stream(write(tmp_path,rows),0,False,150)
    rows[-1]['external_torque_z_mev_mm']=0.
    assert 'orbital_mev_mm' in audit_stream(write(tmp_path,rows),2,False,150)['failed']


def test_prepare_bounded_campaign(tmp_path):
    from af_validate_native_noether import prepare
    m=prepare(tmp_path/'run');assert len(m['cases'])==16
    for c in m['cases']:
        text=Path(c['macro']).read_text()
        assert '/qsim/out/stepPrescale 0' in text
        assert c['first_string_vector_length_mm']<=1e-10
        assert c['nominal_stop_ps']==150
    assert len(list((tmp_path/'run').glob('*.mac')))==16


def test_lorentz_axial_torque_against_off_axis_circle():
    import numpy as np
    # Off-axis orbit has changing r^2; integrate x F_y - y F_x directly.
    theta=np.linspace(.2,1.7,20001);radius=4.;center=np.array([2.,-3.]);k=.299792458*2
    x=center[0]+radius*np.cos(theta);y=center[1]+radius*np.sin(theta)
    dx=-radius*np.sin(theta);dy=radius*np.cos(theta)
    # dP/dtheta = k (dy,-dx) for positive charge in Bz.
    torque=x*(-k*dx)-y*(k*dy)
    integral=np.sum((torque[1:]+torque[:-1])*np.diff(theta)/2)
    exact=-k*((x[-1]**2+y[-1]**2)-(x[0]**2+y[0]**2))/2
    assert abs(integral-exact)<1e-7
