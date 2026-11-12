"""Offline UFO builder tests. Does not validate MadGraph import or matrix elements."""
import ast
from pathlib import Path
import re
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from build_triplet_toy_ufo import EXTRA, REQUIRED_BASE, build_model


class TestTripletBuilder(unittest.TestCase):
    def test_six_distinct_states(self):
        s = EXTRA["particles.py"]
        self.assertEqual(s.count("= Particle("), 3)
        self.assertEqual(s.count(".anti()"), 3)
        self.assertEqual(
            re.findall(r"(?<![a-z])(?:anti)?name='([^']+)'", s),
            ["q1+", "q1-", "q0", "q0~", "q2-", "q2+"],
        )

    def test_eight_vector_vertices(self):
        s = EXTRA["vertices.py"]
        self.assertEqual(s.count("= Vertex("), 8)
        self.assertEqual(s.count("lorentz=[L.FFV1]"), 8)
        self.assertNotIn("P.Q0bar, P.Q0, P.Z", s)

    def test_charge_conservation(self):
        charges = dict(Q1m=-1, Q1p=1, Q0=0, Q0bar=0, Q2m=-1,
                       Q2p=1, a=0, Z=0, W__plus__=1, W__minus__=-1)
        for v in re.findall(r"particles=\[([^]]+)\]", EXTRA["vertices.py"]):
            fields = [part.strip().split(".")[-1] for part in v.split(",")]
            self.assertEqual(sum(charges[p] for p in fields), 0, fields)

    def test_build_and_guard_overwrite(self):
        with tempfile.TemporaryDirectory() as td:
            source = Path(td) / "sm"
            source.mkdir()
            mocks = {
                "particles.py": "import parameters as Param\n# Particle( W__plus__ W__minus__\n",
                "parameters.py": "from object_library import Parameter\n# ee = Parameter, cw = Parameter, sw = Parameter\n",
                "couplings.py": "from object_library import Coupling\n# ee\n",
                "vertices.py": "import particles as P\nimport couplings as C\nimport lorentz as L\n",
                "lorentz.py": 'FFV1 = "Gamma(3,2,1)"\n',
            }
            for f in REQUIRED_BASE:
                (source / f).write_text(mocks.get(f, "\n"))
            dest = Path(td) / "TripletEW_UFO"
            build_model(source, dest)
            for f in REQUIRED_BASE:
                ast.parse((dest / f).read_text())
            with self.assertRaises(FileExistsError):
                build_model(source, dest)
            build_model(source, dest, force=True)


if __name__ == "__main__":
    unittest.main()
