#!/usr/bin/env python3
"""Build NLO-QCD electroweak Dirac triplet UFO from MG5's loop_sm.

The produced UFO is complete for Drell--Yan production at NLO QCD, conditional
on verification with the selected MG5 version. It does NOT contain infracolor
gluons, confinement, weak pion decays or NLO electroweak corrections.

Rationale: triplet is a spectator under SU(3)c. Thus ALL one-loop QCD UV/R2
counterterms belong to the incoming SM-QCD sector, which is imported unchanged
from MG5's validated loop_sm. Do not invent triplet QCD counterterms.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

PROVENANCE_NAME = "QUIRKSIM_NLO_PROVENANCE.json"

ADDITIONS = {
"parameters.py": r'''
# BEGIN QuirkSim Dirac SU(2)_L triplet -- external mass inputs
# Exact tree-level SU(2)_L invariance requires MQC=MQ0. The positive
# charged-neutral O(alpha_EW) mass shift can be applied downstream as a
# pole-mass approximation, but is not part of this NLO-QCD model.
MQ0 = Parameter(name='MQ0',nature='external',type='real',value=225.0,
    texname='m_{Q^0}',lhablock='MASS',lhacode=[9000112])
MQC = Parameter(name='MQC',nature='external',type='real',value=225.0,
    texname='m_{Q^+}',lhablock='MASS',lhacode=[9000111])
# END QuirkSim
''',
"particles.py": r'''
# BEGIN QuirkSim vector-like Dirac triplet (Y=0)
# PDGs are project-local and not assigned by the PDG.
Q1p=Particle(pdg_code=9000111,name='q1+',antiname='q1-',spin=2,
    color=1,mass=Param.MQC,width=Param.ZERO,texname='Q_1^+',
    antitexname='Q_1^-',charge=1,LeptonNumber=0,GhostNumber=0,QuirkNumber=1)
Q1m=Q1p.anti()
Q0=Particle(pdg_code=9000112,name='q0',antiname='q0~',spin=2,
    color=1,mass=Param.MQ0,width=Param.ZERO,texname='Q^0',
    antitexname='\\bar Q^0',charge=0,LeptonNumber=0,GhostNumber=0,QuirkNumber=1)
Q0bar=Q0.anti()
Q2m=Particle(pdg_code=9000113,name='q2-',antiname='q2+',spin=2,
    color=1,mass=Param.MQC,width=Param.ZERO,texname='Q_2^-',
    antitexname='Q_2^+',charge=-1,LeptonNumber=0,GhostNumber=0,QuirkNumber=1)
Q2p=Q2m.anti()
# END QuirkSim
''',
"couplings.py": r'''
# BEGIN QuirkSim vector-like triplet gauge currents
# Vectorlike SU(2)L triplet: T3=diag(+1,0,-1), Y=0.
# W coupling g/sqrt(2)*sqrt(2)=g=ee/sw, both chiralities.
# CT_couplings for NLO QCD are inherited unmodified from loop_sm.
GC_Q1_A=Coupling(name='GC_Q1_A',value='complex(0,1)*ee',
    order={'QED':1})
GC_Q2_A=Coupling(name='GC_Q2_A',value='-complex(0,1)*ee',
    order={'QED':1})
GC_Q1_Z=Coupling(name='GC_Q1_Z',value='complex(0,1)*ee*cw/sw',
    order={'QED':1})
GC_Q2_Z=Coupling(name='GC_Q2_Z',value='-complex(0,1)*ee*cw/sw',
    order={'QED':1})
GC_Q_W=Coupling(name='GC_Q_W',value='complex(0,1)*ee/sw',
    order={'QED':1})
# END QuirkSim
''',
"vertices.py": r'''
# BEGIN QuirkSim: eight FFV vector-like gauge vertices
# The neutral state has no diagonal A/Z coupling at tree-level.
V_Q1_A=Vertex(name='V_Q1_A',particles=[P.Q1m,P.Q1p,P.A],
    color=['1'],lorentz=[L.FFV1],couplings={(0,0):C.GC_Q1_A})
V_Q2_A=Vertex(name='V_Q2_A',particles=[P.Q2p,P.Q2m,P.A],
    color=['1'],lorentz=[L.FFV1],couplings={(0,0):C.GC_Q2_A})
V_Q1_Z=Vertex(name='V_Q1_Z',particles=[P.Q1m,P.Q1p,P.Z],
    color=['1'],lorentz=[L.FFV1],couplings={(0,0):C.GC_Q1_Z})
V_Q2_Z=Vertex(name='V_Q2_Z',particles=[P.Q2p,P.Q2m,P.Z],
    color=['1'],lorentz=[L.FFV1],couplings={(0,0):C.GC_Q2_Z})
V_Q10_Wp=Vertex(name='V_Q10_Wp',particles=[P.Q1m,P.Q0,P.W__plus__],
    color=['1'],lorentz=[L.FFV1],couplings={(0,0):C.GC_Q_W})
V_Q10_Wm=Vertex(name='V_Q10_Wm',particles=[P.Q0bar,P.Q1p,P.W__minus__],
    color=['1'],lorentz=[L.FFV1],couplings={(0,0):C.GC_Q_W})
V_Q02_Wp=Vertex(name='V_Q02_Wp',particles=[P.Q0bar,P.Q2m,P.W__plus__],
    color=['1'],lorentz=[L.FFV1],couplings={(0,0):C.GC_Q_W})
V_Q02_Wm=Vertex(name='V_Q02_Wm',particles=[P.Q2p,P.Q0,P.W__minus__],
    color=['1'],lorentz=[L.FFV1],couplings={(0,0):C.GC_Q_W})
# END QuirkSim
''',
}

def upstream_version(root: Path) -> str|None:
    try:
        return subprocess.check_output(["git","-C",str(root),"rev-parse","HEAD"],
                   text=True,stderr=subprocess.DEVNULL).strip()
    except (OSError,subprocess.CalledProcessError):
        return None

def patch_legacy_param_card_writer(destination: Path) -> bool:
    """Repair Python-2 comparator sorting in old loop_sm parameter-card writers.

    MG5 3.5.1 may still have `need_writing.sort(self.order_param)`, which
    raises TypeError under Python 3. This changes only the *copied* UFO model,
    never the MG5 installation; newer loop_sm writers are left untouched.
    """
    writer = destination / "write_param_card.py"
    source = writer.read_text(encoding="utf-8")
    old = "need_writing.sort(self.order_param)"
    if old not in source:
        return False
    new = "need_writing.sort(key=cmp_to_key(self.order_param))"
    source = source.replace(old, new)
    if "from functools import cmp_to_key" not in source:
        source = "from functools import cmp_to_key\n" + source
    writer.write_text(source, encoding="utf-8")
    return True


def file_sha256(path: Path) -> str:
    """Record upstream file provenance even if MG5 came from a tarball."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(mg5:Path, destination:Path, force:bool=False):
    mg5=mg5.resolve()
    src=mg5/"models"/"loop_sm"
    if not (src/"CT_vertices.py").is_file():
        raise FileNotFoundError(f"Missing NLO SM model: {src}")
    if destination.exists():
        if not force:
            raise FileExistsError(f"{destination} exists; use --force for rebuilding")
        if destination.resolve()==src.resolve():
            raise ValueError("Refusing to overwrite MG5 loop_sm")
        shutil.rmtree(destination)
    shutil.copytree(src,destination)
    compat_patched = patch_legacy_param_card_writer(destination)
    for filename, content in ADDITIONS.items():
        file=destination/filename
        if not file.is_file():
            raise FileNotFoundError(file)
        with file.open("a",encoding="utf-8") as f:
            f.write("\n"+content)
    # Old loop_sm restrict cards are parameter cards for a pure SM model.
    # MG5 automatically imports restrict_default.dat, which must include
    # external MASS entries for *every* BSM mass, otherwise import fails.
    for restriction in destination.glob("restrict_*.dat"):
        contents=restriction.read_text(encoding="utf-8")
        masses=(" 9000111 2.250000e+02 # MQC\n"
                " 9000112 2.250000e+02 # MQ0\n")
        if "Block MASS" not in contents:
            raise ValueError(f"Cannot find MASS block in {restriction}")
        if "9000111" not in contents and "9000112" not in contents:
            contents=contents.replace("Block MASS", "Block MASS\n"+masses,1)
            restriction.write_text(contents,encoding="utf-8")
    # Work in a standard massless-five-flavour initial-state scheme.
    # The old loop_sm default has MB=4.7 GeV; for p including b/bbar
    # NLO massless-parton PDFs this is not a consistent default.
    # Preserve a named 4F/massive-b restriction for explicit studies.
    (destination/"restrict_massive_b.dat").write_text(
        (destination/"restrict_default.dat").read_text(encoding="utf-8"),
        encoding="utf-8")
    (destination/"restrict_default.dat").write_text(
        (destination/"restrict_no_b_mass.dat").read_text(encoding="utf-8"),
        encoding="utf-8")
    manifest={
        "base_model":"MG5_aMC models/loop_sm (upstream NLO QCD CT)",
        "base_git_commit":upstream_version(mg5),
        "base_loop_sm_sha256": {
            filename: file_sha256(src / filename) for filename in
            ("CT_vertices.py", "CT_couplings.py", "CT_parameters.py",
             "vertices.py", "couplings.py", "particles.py",
             "write_param_card.py")
        },
        "legacy_python2_param_card_writer_patched": compat_patched,
        "scope":"pp -> Dirac Y=0 triplet: full LO EW + NLO QCD corrections only",
        "is_fully_dynamical_quirk_model":False,
        "nlo_qcd_uv_r2":"inherited SM QCD counterterms unchanged",
        "initial_flavour_scheme":"massless u d s c b; restrict_no_b_mass is the default",
        "triplet_color":"SU(3)c singlet; hidden SU(N_IC) index not represented",
        "cross_section_hidden_multiplicity":"multiply by N_IC only under spectator hidden-color assumptions",
        "weak_decay_implementation":"not represented; QuirkSim propagator outside UFO",
        "smearing":"charged-neutral splitting must be treated consistently as EW loop-level input",
    }
    (destination/PROVENANCE_NAME).write_text(
        json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print("Built NLO-QCD triplet UFO:",destination)
    print("SM QCD UV/R2 sector inherited from:",src)
    print("Upstream MG5 commit:",manifest["base_git_commit"] or
          "(not a Git checkout; original files identified by SHA256)")
    print("Legacy parameter-card writer patched for Python 3:", compat_patched)
    return destination

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mg5-dir",type=Path,required=True,
                   help="MadGraph5_aMC checkout containing models/loop_sm")
    p.add_argument("--output",type=Path,default=Path("models/DiracTripletEW_NLO_UFO"))
    p.add_argument("--force",action="store_true")
    a=p.parse_args()
    build(a.mg5_dir,a.output,a.force)

if __name__=="__main__":
    main()
