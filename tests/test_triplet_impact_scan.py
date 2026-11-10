"""Guard reproducible nonradial impact scans and initial angular momentum."""
import csv
import json
import math
import subprocess
import sys
from pathlib import Path
import pytest

from quirksim.dynamics import Parameters, HBARC_GEV_CM, propagate
from quirksim.physics import benchmark


def test_nonradial_stable_trajectory_has_nonzero_closest_approach():
    """A finite initial offset and longitudinal momentum imply L != 0."""
    radial=propagate(Parameters(mode='stable',duration_ps=450,impact_cm=0))
    offaxis=propagate(Parameters(mode='stable',duration_ps=450,impact_cm=1e-4))
    assert radial['turning_points']
    assert offaxis['turning_points']
    assert offaxis['turning_points'][0]['separation_cm'] > 1e-7
    assert radial['turning_points'][0]['separation_cm'] < 1e-7
    e=offaxis['segment_energies_gev']
    assert abs(e[-1]-e[0])<1e-4


def test_impact_scan_runner_with_matched_seeds(tmp_path):
    root=Path(__file__).resolve().parents[1]
    cmd=[sys.executable,str(root/'scripts/triplet_impact_scan.py'),
         '--impact-um','0,1','--nevents','2','--duration-ps','450',
         '--seed','42','--progress-every','1','--output-dir',str(tmp_path)]
    result=subprocess.run(cmd,cwd=root,check=True,text=True,
                          capture_output=True,timeout=90)
    assert '[impact-result] b=0 um' in result.stdout
    assert '[impact-result] b=1 um' in result.stdout
    results=json.loads((tmp_path/'summary.json').read_text())
    assert len(results['results'])==2
    assert results['results'][0]['initial_L_over_hbar_abs']==0
    m=benchmark('triplet',225)['charged_mass_gev']
    pz=m*.4/math.sqrt(1-.4*.4)
    assert results['results'][1]['initial_L_over_hbar_abs']==pytest.approx(1e-4*pz/HBARC_GEV_CM)
    with (tmp_path/'paired_events.csv').open() as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==4
    assert set(row['seed'] for row in rows)=={'42','43'}
    assert all(row['channel_matched']=='True' for row in rows)
