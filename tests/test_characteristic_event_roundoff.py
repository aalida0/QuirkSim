import numpy as np
from quirksim.conformal_characteristics import UnitRayHistory
from quirksim.conformal_boundary import massive_boundary


def history():
    h=UnitRayHistory([0.,1.],[[1.,0.,0.]]*2)
    h.jump(1.,[0.,1.,0.])
    h.append(2.,[0.,1.,0.])
    return h


def test_processed_front_and_value_lookup_use_same_roundoff_window():
    h=history()
    for count in (17,24,31):
        z=1.-count*np.spacing(1.)
        assert h.next_jump(z) is None  # event traversal regards this as processed
        np.testing.assert_array_equal(h.value(z),[0.,1.,0.])
        # No zero-step ray defect after storing the post-event reflected ray.
        _,outgoing,_=massive_boundary(np.zeros(3),225.,h.value(z))
        np.testing.assert_array_equal(outgoing,[0.,-1.,0.])


def test_outside_event_window_remains_before_front():
    h=history();z=1.-40*np.spacing(1.)
    assert h.next_jump(z)==1.
    np.testing.assert_array_equal(h.value(z),[1.,0.,0.])
    np.testing.assert_array_equal(h.left_value(1.),[1.,0.,0.])
    np.testing.assert_array_equal(h.value(1.),[0.,1.,0.])
    np.testing.assert_array_equal(h.primitive(1.),[1.,0.,0.])
