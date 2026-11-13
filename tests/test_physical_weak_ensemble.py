"""Inclusive input replay, online minima, censoring and conservation failures."""
import json,subprocess,sys
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from af_physical_weak_ensemble import prepare
from audit_physical_weak import encounter_summary,summarize
from audit_native_noether import audit_stream
from quirksim.weak_plan import make_plan
from test_native_noether import fixture,write


def test_prepare_uses_physical_shared_latents(tmp_path):
    m=prepare(tmp_path/'campaign');assert len(m['cases'])==58
    assert m['independent_events_per_field']==8
    assert len(m['conservation_preflight_cases'])==2
    assert set(m['conservation_preflight_cases']).issubset(c['name'] for c in m['cases'])
    for seed in range(1,9):
        selected=[c for c in m['cases'] if c['seed']==seed]
        expected=make_plan(225.,seed)
        assert all(c['plan']==expected for c in selected)
        assert len({c['general_plan_sha256'] for c in selected})==1
        assert len(selected)==(12 if seed==1 else 6)
        assert all('/qsim/out/stepPrescale 0' in Path(c['macro']).read_text() for c in selected)
    assert make_plan(225.,4)[0]['proper_ps']>980  # Retain long-lifetime censoring; never force it into the window.
    assert m['validation_clocks_ps'] is None


def test_compiled_online_encounter_invariants(tmp_path):
    exe=tmp_path/'encounter'
    subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic',str(ROOT/'research/quirkpheno_recoil/native_encounter_check.cc'),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
    from audit_inclusive_weak import local_minima
    online=json.loads(subprocess.check_output([str(exe),'fixture'],text=True))
    pair=[]
    for i,dt in enumerate((.07,.11)):
        t=np.arange(int(20/dt)+1)*dt
        x=np.column_stack([(.5 if i==0 else -.5)*np.sin(t+(0 if i==0 else .1)),(.03 if i==0 else -.03)*np.cos(2*t),np.zeros(len(t))]);pair.append((t,x))
    offline=local_minima(pair)
    assert len(online)==len(offline)>3
    assert np.allclose(online,[[v['time_cm_ps'],v['separation_mm']] for v in offline],rtol=0,atol=1e-12)



def test_legitimate_nondecay_and_one_decay_are_not_failures(tmp_path):
    # No decay in the sampled finite window is a valid censored history.
    assert audit_stream(write(tmp_path,fixture()),0,True,150,expected_vertices=0)['status'].startswith('PASS')
    rows=fixture();rows.insert(1,dict(rows[0],reason='before_vertex'))
    rows.insert(2,dict(rows[0],reason='after_vertex',vertices=1))
    for d in rows[3:]:d['vertices']=1
    assert audit_stream(write(tmp_path,rows),0,True,150,expected_vertices=1)['vertex_pairs']==1
    rows[-1]['mechanical_p_mev']=[500100,0,0,0]
    assert audit_stream(write(tmp_path,rows),0,True,150,expected_vertices=1)['status']=='FAIL'


def test_censored_encounter_and_invalid_neutral_claim(tmp_path):
    p=tmp_path/'empty.jsonl';p.write_text(json.dumps(dict(kind='end_event',returns=0,processed_intervals=100,buffer_high_water=3))+'\n')
    assert encounter_summary(p,False,0)['status']=='CENSORED_NOT_BOTH_NEUTRAL'
    assert encounter_summary(p,True,100)['status']=='CENSORED_NO_BRACKETED_RETURN'
    e=dict(kind='return',time_cm_ps=105.,separation_mm=.002,endpoints=[dict(pdg=9900017,x_mm=[.001,0,0],p4_mev=[225000,0,0,0]),dict(pdg=-9900017,x_mm=[-.001,0,0],p4_mev=[225000,0,0,0])])
    end=dict(kind='end_event',returns=1,processed_intervals=100,buffer_high_water=3)
    p.write_text(json.dumps(e)+'\n'+json.dumps(end)+'\n')
    assert encounter_summary(p,True,100)['first_return']['separation_mm']==.002
    with pytest.raises(ValueError):encounter_summary(p,False,0)
    with pytest.raises(ValueError):encounter_summary(p,True,110)


def test_cohort_does_not_count_controls_as_independent():
    state={'sampled_seeds':[1],'scope':'test','weak_reference':{}}
    rows=[]
    for field in (0.,2.):
      for mode in ('recoil','collinear_forward','collinear_backward'):
        rows.append(dict(name=f'{field}_{mode}_level8',seed=1,field_tesla=field,level=8,mode=mode,status='PASS_CLASSICAL_PILOT',decay_count=1,decays=[{'channel':'pion'}],geometry={'first_return':None}))
    r=summarize(state,rows)
    assert all(c['independent_events']==1 and c['decay_counts']=={'1':1} for c in r['cohorts'])
    assert all(p['status']=='CENSORED_OR_UNQUALIFIED' for p in r['paired_comparisons'])
    assert r['physical_survival_probability'] is None
