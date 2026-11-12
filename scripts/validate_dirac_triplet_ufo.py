#!/usr/bin/env python3
"""Python-level UFO import and structural validation (NOT MadGraph amplitudes)."""
import sys
import tempfile
from pathlib import Path

def check(folder):
    folder=Path(folder).resolve()
    sys.path.insert(0,str(folder))
    import particles as P
    import parameters as Param
    import vertices as V
    import couplings as C
    import lorentz as L
    from object_library import all_particles, all_vertices
    expected=[
       (P.Q1p,"q1+",9000111,1,1),
       (P.Q1m,"q1-",-9000111,-1,-1),
       (P.Q0,"q0",9000112,0,1),
       (P.Q0bar,"q0~",-9000112,0,-1),
       (P.Q2m,"q2-",9000113,-1,1),
       (P.Q2p,"q2+",-9000113,1,-1),
    ]
    for p,n,pid,charge,qn in expected:
        assert (p.name,p.pdg_code,p.charge,p.QuirkNumber)==(n,pid,charge,qn), p.name
        assert p.spin==2 and p.color==1 and not p.selfconjugate
        assert p.width is Param.ZERO
        assert p in all_particles
    assert len(set(p.pdg_code for p,*_ in expected))==6
    assert Param.MQC.value == Param.MQ0.value == 225.0
    assert L.FFV1.structure == "Gamma(3,2,1)"
    bsm=[V.V_Q1_A,V.V_Q2_A,V.V_Q1_Z,V.V_Q2_Z,
         V.V_Q10_Wp,V.V_Q10_Wm,V.V_Q02_Wp,V.V_Q02_Wm]
    assert len(set(v.name for v in bsm))==8
    for v in bsm:
        assert v in all_vertices and v.lorentz==[L.FFV1] and v.color==['1']
        assert sum(p.charge for p in v.particles)==0,v.name
        assert sum(getattr(p,"QuirkNumber",0) for p in v.particles)==0,v.name
        assert list(v.couplings)==[(0,0)]
    assert [p.name for p in V.V_Q10_Wp.particles]==["q1-","q0","W+"]
    assert [p.name for p in V.V_Q10_Wm.particles]==["q0~","q1+","W-"]
    assert [p.name for p in V.V_Q02_Wp.particles]==["q0~","q2-","W+"]
    assert [p.name for p in V.V_Q02_Wm.particles]==["q2+","q0","W-"]
    assert C.GC_Q1_A.value=="complex(0,1)*ee"
    assert C.GC_Q2_A.value=="-complex(0,1)*ee"
    assert C.GC_Q1_Z.value=="complex(0,1)*ee*cw/sw"
    assert C.GC_Q2_Z.value=="-complex(0,1)*ee*cw/sw"
    assert C.GC_Q_W.value=="complex(0,1)*ee/sw"
    for v in all_vertices:
        if P.Q0 in v.particles or P.Q0bar in v.particles:
            assert P.a not in v.particles and P.Z not in v.particles
    # Actual MG5 3.5.1 bug regression: parameter-card writer must be bundled.
    writer = folder / "write_param_card.py"
    assert writer.is_file(), "Missing write_param_card.py for MG5 3.5.1"
    from write_param_card import ParamCardWriter
    with tempfile.TemporaryDirectory() as td:
        card = Path(td) / "param_card.dat"
        ParamCardWriter(str(card), generic=True)
        text = card.read_text(encoding="utf-8")
        assert "9000111" in text and "9000112" in text and "9000113" in text
        assert "QNUMBERS 9000112" in text, "Neutral Dirac QNUMBERS missing"
    print("PASS: independent parameter card generated with three triplet PDGs")
    print("PASS: 6 Dirac species, 8 EW vertices, 0 neutral-diagonal A/Z")
    print("PASS: all BSM vertices conserve electric charge and quirk number")
    print(f"Loaded {len(all_particles)} particles, {len(all_vertices)} vertices")
    print("NOTE: full MadGraph import, amplitudes, Ward tests still pending")

if __name__=="__main__":
    check(sys.argv[1] if len(sys.argv)>1 else "models/DiracTripletEW_UFO")
