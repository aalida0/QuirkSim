"""Replay only the classical, post-weak-decay kernel from an exact checkpoint."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from quirksim.conformal_characteristics import ConformalString


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkpoint', type=Path)
    p.add_argument('--duration-ps', type=float, default=.01)
    p.add_argument('--max-step-ps', type=float, default=.025)
    p.add_argument('--ray-angle', type=float, default=.0001953125)
    p.add_argument('--causal-fraction', type=float, default=.000390625)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if not np.isfinite(a.duration_ps) or a.duration_ps <= 0:
        p.error('positive finite duration required')
    solver = ConformalString.load_checkpoint(a.checkpoint)
    if np.any(solver.charges != 0):
        p.error('replay requires both endpoints neutral: no future weak-decay scheduling is supplied')
    before = solver.snapshot()
    status = solver.advance(solver.t+a.duration_ps, max_step_ps=a.max_step_ps,
                            max_ray_angle=a.ray_angle, causal_fraction=a.causal_fraction,
                            max_steps=5000, max_rejected_ray_trials=50000,
                            minimum_step_ulps=16, deadline_monotonic=time.monotonic()+30.)
    result = dict(scope='CLASSICAL_POST_DECAY_KERNEL_ONLY_NO_MATCHING_OR_ANNIHILATION',
                  status=status, initial=before, final=solver.snapshot(),
                  diagnostic=solver.last_advance_diagnostic,
                  physical_survival_probability=None)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status=status, diagnostic=solver.last_advance_diagnostic)))
    raise SystemExit(0 if status == 'TARGET_REACHED' else 1)
