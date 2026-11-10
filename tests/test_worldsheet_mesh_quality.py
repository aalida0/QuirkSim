"""Fail closed when the finite-element string collapses / stretches."""
import numpy as np
import pytest
from quirksim.variational_string import StringParameters,VariationalString


def test_initial_equal_spacing_is_acceptable():
    s=VariationalString(StringParameters(nodes=5))
    q=s.mesh_quality(s.initial_state())
    assert q['quality_pass']
    assert q['worst_segment_length_ratio']==pytest.approx(1,abs=1e-12)


def test_degenerate_adjacent_vertices_are_rejected():
    s=VariationalString(StringParameters(nodes=5))
    y=s.initial_state()
    x,v=s._unpack(y)
    x[2]=x[1]
    q=s.mesh_quality(y)
    assert not q['quality_pass']
    assert q['short_relative_to_regulator']


def test_extreme_imbalance_rejected_even_with_nonzero_segment_lengths():
    s=VariationalString(StringParameters(nodes=5))
    y=s.initial_state()
    x,v=s._unpack(y)
    x[1]=x[0]+np.array([1e-6,0,0])
    q=s.mesh_quality(y)
    assert not q['quality_pass']
    assert q['length_ratio_above_50']
