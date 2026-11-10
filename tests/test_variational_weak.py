"""Explicit worldsheet recoil impulse accounting; finite-element gate is essential."""
import math
import numpy as np
import pytest
from quirksim.variational_string import StringParameters,VariationalString


@pytest.mark.parametrize('channel',['pion','electron_neutrino','muon_neutrino'])
def test_weak_decay_canonical_conservation_and_energy_reporting(channel):
    params=StringParameters(nodes=5,initial_beta=.25,initial_length_cm=.05,
                            gauge_mass_fraction=.001)
    s=VariationalString(params)
    before=s.initial_state()
    updated,info=s.weak_decay_impulse(before,0,channel,seed=5,strict=False)
    assert info['channel']==channel
    assert np.linalg.norm(info['momentum_residual_gev'])<1e-7
    assert np.linalg.norm(info['angular_momentum_residual_hbar'])<1e6
    assert math.isfinite(info['energy_residual_gev'])
    assert info['emitted_energy_gev']>.01
    assert s.charge[0]==0 and s.charge[-1]==-1
    assert s.masses[0]==params.mass_gev
    assert not np.array_equal(updated,before)


def test_unacceptably_tight_energy_gate_is_fail_closed():
    s=VariationalString(StringParameters(nodes=5,initial_beta=.3))
    state=s.initial_state()
    initial_mass=float(s.masses[0])
    with pytest.raises(RuntimeError,match='energy mismatch'):
        s.weak_decay_impulse(state,0,'pion',seed=100,energy_tolerance_gev=0)
    assert s.masses[0]==initial_mass
    assert s.charge[0]==1


def test_reject_nonendpoint_or_repeated_decay():
    s=VariationalString(StringParameters())
    state=s.initial_state()
    with pytest.raises(ValueError):
        s.weak_decay_impulse(state,1,'pion')
    final,info=s.weak_decay_impulse(state,0,'pion',seed=4,strict=False)
    with pytest.raises(ValueError):
        s.weak_decay_impulse(final,0,'pion',seed=4,strict=False)
