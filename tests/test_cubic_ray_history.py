import numpy as np
from scipy.integrate import quad_vec
from quirksim.conformal_characteristics import UnitRayHistory, ConformalString
from quirksim.decays import sample_decay


def curve(z):
    q=np.array([1.,.2*np.sin(z),.1*np.cos(2*z)])
    dq=np.array([0.,.2*np.cos(z),-.2*np.sin(2*z)])
    norm=np.linalg.norm(q);n=q/norm
    return n,(dq-n*(n@dq))/norm


def test_cubic_history_converges_against_manufactured_nonplanar_curve():
    errors=[]
    for width in [.4,.2,.1]:
        h=UnitRayHistory([-1.,0.],[curve(-1)[0],curve(0)[0]])
        h.derivatives[-1]=curve(0)[1]
        h.append(width,*curve(width))
        sample=.37*width
        errors.append(np.linalg.norm(h.value(sample)-curve(sample)[0]))
        assert abs(np.linalg.norm(h.value(sample))-1) < 5e-16
        exact=quad_vec(lambda z:curve(z)[0],0.,width,epsabs=1e-14)[0]
        assert np.linalg.norm(h.primitive(width)-h.primitive(0)-exact) < 1e-3*width**4
    assert errors[0]/errors[1]>12 and errors[1]/errors[2]>12


def test_cubic_jump_keeps_distinct_sides_and_has_no_integrated_impulse():
    h=UnitRayHistory([-1.,0.],[[1.,0,0]]*2)
    h.derivatives[-1]=np.zeros(3)
    h.append(1.,[1.,0,0],np.zeros(3));old=h.primitive(1.)
    h.jump(1.,[0.,1.,0.],[.1,0.,0.])
    np.testing.assert_array_equal(h.primitive(1.),old)
    np.testing.assert_array_equal(h.left_value(1.),[1.,0,0])
    np.testing.assert_array_equal(h.value(1.),[0.,1.,0])
    np.testing.assert_allclose(h.gradient(1.,left=True),0.,atol=1e-15)
    np.testing.assert_array_equal(h.gradient(1.),[.1,0,0])


def test_cubic_checkpoint_and_causal_weak_reflection(tmp_path):
    s=ConformalString.boosted_radial(225.2,774.);s.enable_cubic_history()
    for endpoint in (0,1):
        products=sample_decay('pion',s.y[8+3*endpoint:11+3*endpoint],225.2,225.,np.random.default_rng(endpoint))
        s.apply_endpoint_decay(endpoint,225.,products)
    assert s.advance(5.,max_step_ps=.025,causal_fraction=.025)== 'TARGET_REACHED'
    assert s.reflections
    assert s.snapshot()['chord_closure_error_cm']<1e-10
    path=tmp_path/'cubic.npz';s.save_checkpoint(path)
    r=ConformalString.load_checkpoint(path)
    assert r.cubic_history
    for history in ('a','b'):
        a,b=getattr(s,history),getattr(r,history)
        for x,y in zip(a.derivatives,b.derivatives):
            if x is None: assert y is None
            else: np.testing.assert_array_equal(x,y)
    for solver in (s,r): assert solver.advance(5.01,max_step_ps=.005)=='TARGET_REACHED'
    np.testing.assert_array_equal(s.y,r.y)
