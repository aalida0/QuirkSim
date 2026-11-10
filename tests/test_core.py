"""Physics sanity gates; not a substitute for full moving-string validation."""
import math
import numpy as np
import pytest
from quirksim.physics import benchmark, doublet_splitting, triplet_splitting
from quirksim.dynamics import Parameters, pion_two_body, propagate


def test_splittings():
    assert .30 < doublet_splitting(300) < .33
    assert .15 < triplet_splitting(300) < .17


def test_pion_recoil_conservation():
    rng = np.random.default_rng(3)
    for mom in ([0,0,0],[200,0,75],[0,0,-300]):
        mass_parent, mass_daughter = 225.16,225.
        daughter,pion=pion_two_body(np.array(mom,float),mass_parent,mass_daughter,rng)
        assert np.allclose(daughter+pion[1],mom,atol=1e-9)
        assert math.isclose(math.hypot(mass_parent,np.linalg.norm(mom)),
                            math.hypot(mass_daughter,np.linalg.norm(daughter))+pion[0],
                            abs_tol=1e-9)


def test_triplet_decay():
    b=benchmark("triplet",300)
    assert 100 < b["lifetime_ps"] < 500
    assert .9 < b["pion_branching"] < 1


def test_stable_conservation():
    result=propagate(Parameters(mode="stable",duration_ps=40,impact_cm=.002))
    energies=np.array(result["segment_energies_gev"])
    assert np.max(np.abs(energies-energies[0])) < 5e-5
    assert not result["events"]


def test_decay_propagation():
    result=propagate(Parameters(mode="recoil",duration_ps=700,seed=8))
    assert len(result["events"]) <= 2
    assert result["end_time_ps"] <= 700
    assert result["maximum_separation_cm"] > result["minimum_separation_including_production_cm"]


def test_input_validation():
    with pytest.raises(ValueError):
        Parameters(lambda_ev=0)
    with pytest.raises(ValueError):
        Parameters(mode="unknown")
