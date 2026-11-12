"""Tests for reproducible fixed-order 6-channel triplet NLO card generation."""
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from prepare_triplet_nlo_inclusive import make_cards


class TestTripletNLOInclusiveCards(unittest.TestCase):
    def fixture(self, root):
        model = root / "DiracTripletEW_NLO_UFO"
        model.mkdir()
        (model / "CT_vertices.py").write_text("# NLO QCD counterterms\n")
        return model

    def test_internal_pdf_pilot(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            model=self.fixture(root)
            generate, launch=make_cards(
                model=model,output_dir=root/"nlo",card_dir=root/"cards",
            )
            a=generate.read_text()
            b=launch.read_text()
            self.assertEqual(a.count(" QCD=0 QED=2 [QCD] @"),6)
            self.assertEqual(a.count("generate p p >"),1)
            self.assertEqual(a.count("add process p p >"),5)
            self.assertIn("define p = g u u~ d d~ s s~ c c~ b b~",a)
            self.assertTrue(a.index("generate")<a.index("output"))
            self.assertIn("fixed_order=on",b)
            self.assertIn("set req_acc_FO 0.02",b)
            self.assertIn("set ebeam1 6800",b)
            self.assertIn("set ebeam2 6800",b)
            self.assertNotIn("set pdlabel",b)
            self.assertNotIn("set lhaid",b)
            (root/"nlo").mkdir()
            with self.assertRaises(FileExistsError):
                make_cards(model=model,output_dir=root/"nlo",
                           card_dir=root/"other_cards")
    def test_lhapdf_mode_card(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            model=self.fixture(root)
            _, launch=make_cards(
                model=model,output_dir=root/"nlo",card_dir=root/"cards",
                pdf="lhapdf",lhaid=303400,accuracy=0.01,
            )
            self.assertIn("set pdlabel lhapdf",launch.read_text())
            self.assertIn("set lhaid 303400",launch.read_text())
            self.assertIn("set req_acc_FO 0.01",launch.read_text())
    def test_invalid_modes(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            model=self.fixture(root)
            for opts in ({"pdf":"lhapdf"},{"pdf":"internal","lhaid":303400},
                         {"accuracy":0},{"accuracy":float("nan")},
                         {"lhaid":0,"pdf":"lhapdf"}):
                with self.subTest(opts=opts):
                    with self.assertRaises(ValueError):
                        make_cards(model=model,output_dir=root/"nlo",
                                   card_dir=root/"cards",**opts)
            with self.assertRaises(FileNotFoundError):
                make_cards(model=root/"missing",output_dir=root/"nlo",
                           card_dir=root/"cards")


if __name__=="__main__":
    unittest.main()
