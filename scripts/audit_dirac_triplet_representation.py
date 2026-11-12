#!/usr/bin/env python3
"""Independent SU(2)L representation audit of the Dirac triplet UFO currents."""
import cmath
import math
import numpy as np


def check():
    # Isospin-1 representation; q1+, q0, q2- states with T3=+1,0,-1.
    t3=np.diag([1.,0.,-1.]).astype(complex)
    raising=np.array([[0.,math.sqrt(2),0.],
                      [0.,0.,math.sqrt(2)],[0.,0.,0.]],dtype=complex)
    lowering=raising.conj().T
    assert np.allclose(raising@lowering-lowering@raising,2*t3)
    assert np.allclose(t3@raising-raising@t3,raising)
    assert np.allclose(t3@lowering-lowering@t3,-lowering)
    q=t3.copy()  # Y=0.
    sw=0.48
    cw=math.sqrt(1-sw**2)
    e=0.31
    gw=e/sw
    gA=e*q
    gZ=gw*cw*t3
    gWp=gw/math.sqrt(2)*raising
    assert np.allclose(np.diag(gA),[e,0.,-e])
    assert np.allclose(np.diag(gZ),[e*cw/sw,0.,-e*cw/sw])
    assert abs(gWp[0,1]-e/sw)<1e-14
    assert abs(gWp[1,2]-e/sw)<1e-14
    assert not np.any(gWp-np.triu(gWp))
    print("PASS: SU(2) triplet generator algebra, all EW current coefficients, Y=0")


if __name__=="__main__":
    check()
