import json
import numpy as np
from quirksim.conformal_characteristics import ConformalString
from quirksim.conformal_evolution import evolve


def test_budget_stop_records_state_without_changing_accepted_trajectory():
    limited = ConformalString.radial(225., 774.)
    reference = ConformalString.radial(225., 774.)
    assert limited.advance(.1, max_step_ps=.01, max_steps=1) == 'STOP_STEP_BUDGET'
    assert reference.advance(.01, max_step_ps=.01) == 'TARGET_REACHED'
    np.testing.assert_array_equal(limited.y, reference.y)
    assert limited.a.z == reference.a.z and limited.b.z == reference.b.z
    d = limited.last_advance_diagnostic
    assert d['time_ps'] == .01 and d['target_ps'] == .1
    assert d['attempted_steps'] == 2 and d['max_steps'] == 1
    json.dumps(d, allow_nan=False)
    assert limited.advance(.02, max_step_ps=.01) == 'TARGET_REACHED'
    assert limited.last_advance_diagnostic is None


def test_failed_evolution_exposes_diagnostic_without_quantum_survival():
    r = evolve(duration_ps=.1, max_steps_per_advance=1, max_step_ps=.01)
    assert r['termination'] == 'STOP_STEP_BUDGET'
    assert r['last_advance_diagnostic']['reason'] == r['termination']
    assert r['physical_survival_probability'] is None
    json.dumps(r, allow_nan=False)


def test_optional_sampled_closure_stop_remains_a_numerical_failure():
    r = evolve(duration_ps=1., closure_stop_in_reference_units=1e-20)
    assert r['termination'] == 'STOP_SAMPLED_CLOSURE_TOLERANCE'
    assert r['last_advance_diagnostic']['time_ps'] == r['final']['time_ps']
    assert r['physical_survival_probability'] is None
