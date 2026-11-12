#!/usr/bin/env python3
"""Check spectator Dirac triplet NLO QCD model and inherited SM counterterms."""
from __future__ import annotations
import argparse
import importlib
import math
from pathlib import Path
import sys
import tempfile


def check(mg5_dir:Path, model_dir:Path):
    mg5_dir=mg5_dir.resolve()
    model_dir=model_dir.resolve()
    if not (model_dir/"CT_vertices.py").is_file():
        raise ValueError("Missing inherited NLO QCD UV/R2 terms")
    sys.path.insert(0,str(mg5_dir/"models"))
    sys.path.insert(0,str(model_dir.parent))
    base=importlib.import_module("loop_sm")
    mod=importlib.import_module(model_dir.name)
    P,C,V,Par,L=mod.particles,mod.couplings,mod.vertices,mod.parameters,mod.lorentz
    from pathlib import Path
    expected=[
        (P.Q1p,9000111,1,1),(P.Q1m,-9000111,-1,-1),
        (P.Q0,9000112,0,1),(P.Q0bar,-9000112,0,-1),
        (P.Q2m,9000113,-1,1),(P.Q2p,-9000113,1,-1),
    ]
    for obj,pdg,charge,qn in expected:
        assert (obj.pdg_code,obj.charge,obj.QuirkNumber)==(pdg,charge,qn)
        assert obj.color==1 and obj.spin==2 and obj.width is Par.ZERO
    assert Par.MQ0.value==Par.MQC.value==225.0
    assert P.Q1p.mass is Par.MQC and P.Q2m.mass is Par.MQC
    assert P.Q0.mass is Par.MQ0
    assert len(mod.all_particles)==len(base.all_particles)+6
    assert len(mod.all_vertices)==len(base.all_vertices)+8
    assert len(mod.all_CTvertices)==len(base.all_CTvertices)
    # The UFO module exports one shared all_couplings registry for both
    # tree and CT couplings; five new EW tree couplings are expected.
    assert len(mod.all_couplings)==len(base.all_couplings)+5
    assert len(mod.CT_parameters.all_CTparameters)==len(base.CT_parameters.all_CTparameters)
    assert any("UV" in v.name or "R2" in v.name for v in mod.all_CTvertices)
    for bsm in (V.V_Q1_A,V.V_Q2_A,V.V_Q1_Z,V.V_Q2_Z,
                V.V_Q10_Wp,V.V_Q10_Wm,V.V_Q02_Wp,V.V_Q02_Wm):
        assert bsm in mod.all_vertices
        assert bsm.color==['1']
        assert bsm.lorentz==[L.FFV1] and L.FFV1.structure=='Gamma(3,2,1)'
        assert sum(p.charge for p in bsm.particles)==0
        assert sum(getattr(p,'QuirkNumber',0) for p in bsm.particles)==0
        assert all(getattr(p,'color',1)==1 for p in bsm.particles if p.spin==2)
    assert all(getattr(p,"pdg_code",0) not in [9000111,9000112,9000113]
               for v in mod.all_CTvertices for p in v.particles)
    assert [x.name for x in V.V_Q10_Wp.particles]==['q1-','q0','W+']
    assert [x.name for x in V.V_Q10_Wm.particles]==['q0~','q1+','W-']
    assert [x.name for x in V.V_Q02_Wp.particles]==['q0~','q2-','W+']
    assert [x.name for x in V.V_Q02_Wm.particles]==['q2+','q0','W-']
    assert all(P.Q0 not in v.particles and P.Q0bar not in v.particles
               for v in [V.V_Q1_A,V.V_Q2_A,V.V_Q1_Z,V.V_Q2_Z])
    assert L.FFV1.structure=='Gamma(3,2,1)'
    couplings=[
        (C.GC_Q1_A,1j),(C.GC_Q2_A,-1j),
        (C.GC_Q1_Z,1j*0.88/0.47),(C.GC_Q2_Z,-1j*0.88/0.47),
        (C.GC_Q_W,1j/0.47)
    ]
    context={'complex':complex,'ee':1.,'sw':.47,'cw':.88}
    for obj,expected_value in couplings:
        actual=eval(obj.value,{'__builtins__':{}},context)
        assert abs(actual-expected_value)<1e-12,(obj.name,actual,expected_value)
        assert obj.order=={'QED':1}
    assert mod.coupling_orders.QCD.perturbative_expansion==1
    assert mod.coupling_orders.QED.perturbative_expansion==0
    with tempfile.TemporaryDirectory() as tmp:
        card=Path(tmp)/"param_card.dat"
        # loop_sm ships a legacy standalone writer with absolute imports;
        # supply the actual packaged parameter objects explicitly.
        ext=[p for p in mod.all_parameters if p.nature=='external']
        mod.write_param_card.ParamCardWriter(str(card),list_of_parameters=ext,generic=False)
        card_text=card.read_text()
        for code in (9000111,9000112):
            assert str(code) in card_text
    print(f"PASS: 6 triplet/anti-triplet species, 8 EW vector vertices, vectorlike W coupling")
    print(f"PASS: SM NLO QCD CT vertices retained ({len(mod.all_CTvertices)})")
    print("PASS: NLO QCD coupling order and generated parameter card")
    print("NOTE: NLO process-generation and integrated cross-section validation are separate")

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--mg5-dir",required=True,type=Path)
    parser.add_argument("--model",type=Path,default=Path("models/DiracTripletEW_NLO_UFO"))
    args=parser.parse_args()
    check(args.mg5_dir,args.model)
