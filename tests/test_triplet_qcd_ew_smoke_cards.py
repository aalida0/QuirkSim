"""Ensure experimental mixed NLO QCD+EW generator cannot run with QCD-only UFO."""
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from prepare_triplet_qcd_ew_smoke import prepare
from audit_triplet_nlo_ew_candidate import REQUIRED_FILES

class TestEWGenerationGate(unittest.TestCase):
    def make_ufo(self,root,include_ew):
        p=root/"model"
        p.mkdir()
        for x in REQUIRED_FILES:
            (p/x).write_text("# test\n")
        (p/"particles.py").write_text("Qp1=object()\nQzero=object()\n")
        text="CTVertex(name='uv',type='UV',particles=[P.Qp1,P.Qp1,P.A])\n"
        if include_ew:
            text+="CTVertex(name='r2',type='R2',particles=[P.Qp1,P.Qzero,P.W])\n"
        (p/"CT_vertices.py").write_text(text)
        return p

    def test_reject_without_acknowledge(self):
        with tempfile.TemporaryDirectory() as d:
            base=Path(d)
            p=self.make_ufo(base,True)
            with self.assertRaisesRegex(ValueError,"Acknowledge"):
                prepare(p,base/"out",base/"cards")

    def test_reject_qcd_only(self):
        with tempfile.TemporaryDirectory() as d:
            base=Path(d)
            p=self.make_ufo(base,False)
            with self.assertRaisesRegex(ValueError,"lacks triplet-specific"):
                prepare(p,base/"out",base/"cards",acknowledge_experimental=True)

    def test_experimental_cards_include_photon_and_mixed_corrections(self):
        with tempfile.TemporaryDirectory() as d:
            base=Path(d)
            p=self.make_ufo(base,True)
            g,l=prepare(p,base/"out",base/"cards",
                        acknowledge_experimental=True)
            source=g.read_text()
            self.assertIn("[QCD QED]",source)
            self.assertIn("q1+ q1-",source)
            self.assertIn("q2- q2+",source)
            self.assertIn(" b b~ a",source)
            self.assertIn("set complex_mass_scheme True",source)
            self.assertIn("set lhaid 335900",l.read_text())


if __name__=="__main__":
    unittest.main()
