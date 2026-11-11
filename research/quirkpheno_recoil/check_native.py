"""Compile the independent C++ core and compare to QuirkSim's physical pion sampler."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from quirksim.decays import sample_decay
from quirksim.physics import benchmark
class FixedDirection:
    def __init__(self,sign):self.count=0;self.sign=sign
    def uniform(self,a,b):
        self.count+=1
        return 0 if self.count==1 else (np.pi/2 if self.sign>0 else 3*np.pi/2)
def check(output=None):
    here=Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory(prefix='qp-native-') as tmp:
        exe=Path(tmp)/'check'
        subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic',str(here/'native_check.cc'),'-o',str(exe)],check=True)
        fixtures=json.loads(subprocess.check_output([str(exe)],text=True))
        clock=Path(tmp)/'clock'
        subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic',str(here/'native_clock_check.cc'),'-o',str(clock)],check=True)
        subprocess.run([str(clock)],check=True)
    assert abs(benchmark('triplet',225.)['charged_mass_gev']-225.1571702763212)<1e-12
    maximum=0.
    for fixture in fixtures:
        py=sample_decay('pion',np.array(fixture['parent'][1:])/1000,225.1571702763212,225.,FixedDirection(fixture['sign']))
        for key in ('daughter','pion'):
            maximum=max(maximum,float(max(abs(py[key]-np.array(fixture[key])/1000))))
    assert maximum<3e-9
    if output:Path(output).write_text(json.dumps(fixtures,indent=2)+'\n')
    print('PASS native handoff/core, eight independent Python comparisons; maximum difference GeV',maximum)
    return fixtures
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output');a=parser.parse_args();check(a.output)
