"""Smoke-test the streaming ensemble runner and its checkpoint CSV format."""
import csv
import json
import subprocess
import sys
from pathlib import Path


def test_progress_checkpoint_and_matched_modes(tmp_path):
    root=Path(__file__).resolve().parents[1]
    cmd=[sys.executable,str(root/'scripts/paired_compare.py'),
         '--model','triplet','--mass','225','--lambda-ev','774',
         '--nevents','2','--duration-ps','300','--seed','100',
         '--progress-every','1','--checkpoint-every','1',
         '--skip-stable','--output-dir',str(tmp_path)]
    result=subprocess.run(cmd,cwd=root,text=True,capture_output=True,check=True,timeout=90)
    assert '[QuirkSim] 1/2 matched events' in result.stdout
    with (tmp_path/'events.csv').open(newline='') as f:
        rows=list(csv.DictReader(f))
    with (tmp_path/'pairs.csv').open(newline='') as f:
        pairs=list(csv.DictReader(f))
    assert len(rows)==4
    assert len(pairs)==2
    assert {r['mode'] for r in rows}=={'no_recoil','recoil'}
    summary=json.loads((tmp_path/'summary.json').read_text())
    assert summary['benchmark']['n']==2
    assert set(summary['mode_summaries'])=={'no_recoil','recoil'}
