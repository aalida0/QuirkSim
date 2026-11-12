#!/usr/bin/env python3
"""Create an LO electroweak-only Dirac Y=0 triplet UFO by extending MG5's SM UFO.

No infracolor force, confinement, weak decays, showers, or detector effects.
Run against the *installed MG5* SM UFO to preserve its complete SM definitions.
"""
from __future__ import annotations

import argparse
import ast
import os
from pathlib import Path
import shutil

MODEL_NAME = "DiracTripletEW_UFO"
REQUIRED_BASE = (
    "__init__.py", "object_library.py", "particles.py", "parameters.py",
    "couplings.py", "vertices.py", "lorentz.py", "coupling_orders.py",
    "function_library.py",
)

# '+' and '-' form part of MadGraph particle names; q1- != q2-.
EXTRA = {
    "parameters.py": r'''
# --- QuirkSim: Dirac (SU(2)_L,Y)=(3,0) electroweak-only benchmark ---
# The small EW one-loop splitting is given as an INPUT, not calculated here.
MQ0 = Parameter(name='MQ0', nature='external', type='real',
    value=225.0, texname='m_{Q^0}', lhablock='MASS', lhacode=[9000112])
MQC = Parameter(name='MQC', nature='external', type='real',
    value=225.0, texname='m_{Q^\\pm}', lhablock='MASS', lhacode=[9000111])
''',
    "particles.py": r'''
# --- QuirkSim charged/neutral Dirac triplet (color=1 is a PROXY) ---
# PDG assignments are custom to this toy and NOT official PDG registrations.
Q1p = Particle(pdg_code=9000111, name='q1+', antiname='q1-', spin=2,
    color=1, mass=Param.MQC, width=Param.ZERO, texname='Q_1^+',
    antitexname='Q_1^-', charge=1, GhostNumber=0, LeptonNumber=0, QuirkNumber=1, Y=0)
Q1m = Q1p.anti()
Q0 = Particle(pdg_code=9000112, name='q0', antiname='q0~', spin=2,
    color=1, mass=Param.MQ0, width=Param.ZERO, texname='Q^0',
    antitexname='\\bar{Q}^0', charge=0, GhostNumber=0, LeptonNumber=0, QuirkNumber=1, Y=0)
Q0bar = Q0.anti()
Q2m = Particle(pdg_code=9000113, name='q2-', antiname='q2+', spin=2,
    color=1, mass=Param.MQC, width=Param.ZERO, texname='Q_2^-',
    antitexname='Q_2^+', charge=-1, GhostNumber=0, LeptonNumber=0, QuirkNumber=1, Y=0)
Q2p = Q2m.anti()
''',
    "couplings.py": r'''
# --- QuirkSim vectorlike triplet gauge couplings (both chiralities) ---
# L_int = e A_mu( Q1bar gamma Q1 - Q2bar gamma Q2)
#       + g*cW Z_mu( Q1bar gamma Q1 - Q2bar gamma Q2)
#       + g W+_mu( Q1bar gamma Q0 + Q0bar gamma Q2 ) + h.c.
# g=ee/sw.  The SU(2) triplet Clebsch sqrt(2) cancels W sqrt(2).
GC_Q1_A = Coupling(name='GC_Q1_A', value='complex(0,1)*ee',
    order={'QED':1})
GC_Q2_A = Coupling(name='GC_Q2_A', value='-complex(0,1)*ee',
    order={'QED':1})
GC_Q1_Z = Coupling(name='GC_Q1_Z', value='complex(0,1)*ee*cw/sw',
    order={'QED':1})
GC_Q2_Z = Coupling(name='GC_Q2_Z', value='-complex(0,1)*ee*cw/sw',
    order={'QED':1})
GC_Q_W = Coupling(name='GC_Q_W', value='complex(0,1)*ee/sw',
    order={'QED':1})
''',
    "vertices.py": r'''
# --- QuirkSim tree-level vector gauge vertices ---
# FFV1 is Gamma(3,2,1), not a purely left-handed projector.
V_Q1_A = Vertex(name='V_Q1_A', particles=[P.Q1m, P.Q1p, P.a],
    color=['1'], lorentz=[L.FFV1], couplings={(0,0):C.GC_Q1_A})
V_Q2_A = Vertex(name='V_Q2_A', particles=[P.Q2p, P.Q2m, P.a],
    color=['1'], lorentz=[L.FFV1], couplings={(0,0):C.GC_Q2_A})
V_Q1_Z = Vertex(name='V_Q1_Z', particles=[P.Q1m, P.Q1p, P.Z],
    color=['1'], lorentz=[L.FFV1], couplings={(0,0):C.GC_Q1_Z})
V_Q2_Z = Vertex(name='V_Q2_Z', particles=[P.Q2p, P.Q2m, P.Z],
    color=['1'], lorentz=[L.FFV1], couplings={(0,0):C.GC_Q2_Z})
V_Q10_Wp = Vertex(name='V_Q10_Wp', particles=[P.Q1m, P.Q0, P.W__plus__],
    color=['1'], lorentz=[L.FFV1], couplings={(0,0):C.GC_Q_W})
V_Q10_Wm = Vertex(name='V_Q10_Wm', particles=[P.Q0bar, P.Q1p, P.W__minus__],
    color=['1'], lorentz=[L.FFV1], couplings={(0,0):C.GC_Q_W})
V_Q02_Wp = Vertex(name='V_Q02_Wp', particles=[P.Q0bar, P.Q2m, P.W__plus__],
    color=['1'], lorentz=[L.FFV1], couplings={(0,0):C.GC_Q_W})
V_Q02_Wm = Vertex(name='V_Q02_Wm', particles=[P.Q2p, P.Q0, P.W__minus__],
    color=['1'], lorentz=[L.FFV1], couplings={(0,0):C.GC_Q_W})
''',
}


def build_model(sm_ufo: Path, out: Path, *, force: bool = False) -> Path:
    """Build a self-contained model by copying MG5's SM UFO then appending BSM states.

    The source model is not edited. This is a tree-level illustrative tool.
    """
    sm_ufo = sm_ufo.resolve()
    out = out.resolve()
    if sm_ufo == out or sm_ufo in out.parents:
        raise ValueError('Output cannot be inside the source UFO')
    missing = [p for p in REQUIRED_BASE if not (sm_ufo / p).is_file()]
    if missing:
        raise FileNotFoundError(f"Not an SM UFO; missing {missing} in {sm_ufo}")
    if out.exists():
        if not force:
            raise FileExistsError(f"Output exists: {out}; use --force to replace")
        if not (out / '.quirksim_toy_ufo_marker').is_file():
            raise ValueError('Refusing --force on a directory not built by this script')
        shutil.rmtree(out)
    out.mkdir(parents=True)
    # Copy only UFO code and the standard SM parameter card, not generated caches,
    # decay tables, SM restriction cards, or MG5-specific restrict helper files.
    for f in sm_ufo.iterdir():
        if f.suffix == '.py' and f.name not in ('decays.py', 'build_restrict.py', 'write_param_card.py'):
            shutil.copy2(f, out / f.name)
        elif f.name == 'param_card.dat':
            # Card is kept only as a reference. MadGraph regenerates the authoritative
            # card from external Parameter definitions when importing the model.
            shutil.copy2(f, out / 'sm_reference_param_card.dat')
    for f, patch in EXTRA.items():
        target = out / f
        source = target.read_text(encoding='utf-8')
        # Verify this particular SM UFO has the same expected conventions.
        expected = {
            'particles.py': ('import parameters as Param', 'Particle('),
            'parameters.py': ('from object_library', 'ee = Parameter', 'cw = Parameter', 'sw = Parameter'),
            'couplings.py': ('Coupling', 'ee'),
            'vertices.py': ('import particles as P', 'import couplings as C', 'import lorentz as L'),
        }[f]
        if not all(t in source for t in expected):
            raise ValueError(f'Unsupported SM UFO conventions in {f}')
        target.write_text(source + '\n' + patch + '\n', encoding='utf-8')
        ast.parse(target.read_text(encoding='utf-8'), filename=f)
    # Verify essential Lorentz definitions exist and anti-particle fields are available.
    for f, tokens in {
        'lorentz.py': ('FFV1', 'Gamma(3,2,1)'),
        'particles.py': ('W__plus__', 'W__minus__', "name='q1+'", "name='q2-'"),
    }.items():
        s=(out/f).read_text(encoding='utf-8')
        if not all(token in s for token in tokens):
            raise ValueError(f'Unsupported SM UFO conventions: missing {tokens} in {f}')
    (out / '.quirksim_toy_ufo_marker').write_text('Generated by QuirkSim toy UFO builder.\n')
    (out / 'README_QUIRKSIM.txt').write_text(
        'QuirkSim TripletEW_UFO: LO EW-only Dirac triplet Y=0; color-singlet proxy.\n'
        'States q1+/q1-, q0/q0~, q2-/q2+. External masses MQC, MQ0.\n'
        'NO infracolor, strings, radiation, weak pion decays or detector.\n'
        'Do not interpret the generated cross section as a physical IC prediction.\n'
        'Reference generator: SM UFO from installed MadGraph5_aMC.\n'
    )
    return out


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sm-ufo', type=Path, required=True,
                   help='SM UFO directory (usually MG5_aMC*/models/sm)')
    p.add_argument('--out', type=Path, required=True, help='Output UFO directory')
    p.add_argument('--force', action='store_true', help='Replace only previously generated output')
    a=p.parse_args()
    dest=build_model(a.sm_ufo,a.out,force=a.force)
    print('Built:', dest)
    print('Now run: import model ' + dest.name + ' (from MadGraph models/ directory)')


if __name__ == '__main__':
    main()
