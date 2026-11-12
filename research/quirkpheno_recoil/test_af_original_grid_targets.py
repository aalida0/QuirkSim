"""Prepare-only regression checks; no Geant4 installation or job submission."""
import hashlib
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from af_original_grid_targets import prepare, TARGETS


def test_default_pilot_cards_and_paired_plans(tmp_path):
    state = prepare(tmp_path / 'pilot')
    assert len(state['cases']) == 102
    assert state['numerical_validation'] == 'NOT_RUN'
    assert state['full_weak_channels_supported'] is False
    assert {(x['mass_gev'], x['lambda_ev']) for x in state['points']} == set(TARGETS)
    plans = {}
    for case in state['cases']:
        macro = Path(case['macro']).read_text()
        assert f"/quirk/mass {case['mass_gev']:.17g} GeV" in macro
        assert f"/quirk/lambda {case['lambda_ev']:.17g} eV" in macro
        assert case['mass_gev'] > case['neutral_mass_gev'] + .13957039
        assert hashlib.sha256(macro.encode()).hexdigest() == case['macro_sha256']
        assert '/qsim/maxTime 4 ns' in macro
        if not case['weak_enabled']:
            continue
        key = case['point'], case['seed']
        if key in plans:
            assert plans[key] == case['plan']
        plans[key] = case['plan']
        raw = Path(case['plan_file']).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == case['plan_sha256']
        for endpoint in case['plan']:
            assert endpoint['proper_ps'] > 0
            assert math.isclose(np.linalg.norm(endpoint['direction']), 1., abs_tol=1e-12)
    for point in state['points']:
        assert 0 < point['pion_pair_fraction'] <= 1
        assert math.isclose(point['pion_pair_fraction'] + point['omitted_channel_pair_fraction'], 1.)

