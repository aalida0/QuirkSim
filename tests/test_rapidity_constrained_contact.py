"""Regression: bounded-speed constrained NG coordinates and runnable CLI."""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from quirksim.equidistributed_ng import EquidistributedNGString
from quirksim.variational_string import StringParameters
from quirksim.variational_rapidity import physical_to_rapidity,rapidity_to_physical


def run_short_cli(tmp_path,nodes=5,duration=4.,regulator=1e-9,step=.1):
    root=Path(__file__).resolve().parents[1]
    out=tmp_path/f'short_{nodes}_{duration}.json'
    cmd=[sys.executable,str(root/'scripts/rapidity_constrained_contact.py'),
         '--nodes',str(nodes),'--duration-ps',str(duration),
         '--regulator-cm',str(regulator),'--pinch-fraction','0.0004',
         '--max-step-ps',str(step),'--output',str(out)]
    subprocess.run(cmd,cwd=root,check=True,text=True,
                   capture_output=True,timeout=90)
    return json.loads(out.read_text())


def test_rapidity_constrained_coordinate_roundtrip_for_zero_mass_gauge():
    s=EquidistributedNGString(StringParameters(nodes=5,
                           initial_beta=.4,initial_length_cm=.1))
    state=s.boosted_straight_initial_state(beta_cm_x=.25)
    w=physical_to_rapidity(state,s.n)
    recovered=rapidity_to_physical(w,s.n)
    assert np.max(np.abs(recovered-state))<1e-12


def test_short_field_trajectory_matches_original_coordinate(tmp_path):
    params=StringParameters(nodes=5,initial_beta=.4,initial_length_cm=.1,
                            b_tesla=2.,regularizer_cm=1e-9,max_step_ps=.1)
    s=EquidistributedNGString(params)
    y=s.boosted_straight_initial_state(beta_cm_x=.25)
    base,sol=s.solve(y,duration_ps=4.)
    x,v=s._unpack(sol.y[:,-1])
    rapid=run_short_cli(tmp_path)
    assert rapid['status']=='REACHED_HORIZON'
    assert rapid['max_vertex_beta']<1
    assert rapid['max_gauge_length_deviation_cm']<1e-7
    assert rapid['energy_drift_gev']==pytest.approx(base['energy_drift_gev'],abs=1e-7)
    assert abs(rapid['endpoint_chord_cm']-
               float(np.linalg.norm(x[-1]-x[0])))<1e-7


def test_short_run_refuses_unsupported_physical_claims(tmp_path):
    result=run_short_cli(tmp_path,nodes=3,duration=.05,step=.05)
    assert result['status']=='REACHED_HORIZON'
    assert 'NOT validated continuum NG' in result['physics_scope']
