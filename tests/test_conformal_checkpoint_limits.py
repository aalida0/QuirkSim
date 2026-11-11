import time
import numpy as np
from quirksim.conformal_characteristics import ConformalString
from quirksim.decays import sample_decay


def test_checkpoint_preserves_jumps_integrals_and_continuation(tmp_path):
    s = ConformalString.boosted_radial(225.2, 774.)
    for endpoint in (0, 1):
        products = sample_decay('pion', s.y[8+3*endpoint:11+3*endpoint],
                                225.2, 225., np.random.default_rng(endpoint))
        s.apply_endpoint_decay(endpoint, 225., products)
    assert s.advance(.02, max_step_ps=.01) == 'TARGET_REACHED'
    file = tmp_path/'checkpoint.npz'
    s.save_checkpoint(file)
    restored = ConformalString.load_checkpoint(file)
    for name in ('a', 'b'):
        for field in ('z', 'n', 'integrals', 'jumps'):
            np.testing.assert_array_equal(getattr(getattr(s, name), field),
                                          getattr(getattr(restored, name), field))
    assert restored.emissions == s.emissions
    for solver in (s, restored):
        assert solver.advance(.04, max_step_ps=.01) == 'TARGET_REACHED'
    np.testing.assert_array_equal(s.y, restored.y)
    assert s.snapshot() == restored.snapshot()


def test_deadline_and_resolution_guards_leave_unaccepted_state_unchanged():
    for controls, reason in [
        ({'deadline_monotonic': time.monotonic()-1}, 'STOP_DIAGNOSTIC_WALL_TIME'),
        ({'max_step_ps': 1e-17, 'minimum_step_ulps': 16}, 'STOP_DIAGNOSTIC_TIME_RESOLUTION')]:
        s = ConformalString.radial(225., 774.)
        initial = s.y.copy()
        assert s.advance(1., **controls) == reason
        np.testing.assert_array_equal(s.y, initial)
        assert s.t == 0 and s.last_advance_diagnostic['reason'] == reason


def test_rejected_ray_budget_is_a_failure_not_accepted_step():
    s = ConformalString.boosted_radial(225., 774., b_tesla=1000.)
    initial = s.y.copy()
    assert s.advance(.1, max_ray_angle=1e-10,
                     max_rejected_ray_trials=1) == 'STOP_DIAGNOSTIC_REJECTION_BUDGET'
    np.testing.assert_array_equal(s.y, initial)
    assert s.last_advance_diagnostic['rejected_ray_trials'] == 1
