"""Synthetic streaming checks for large NLOCT UV/R2 content inspection."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from inspect_triplet_nloct_large import inspect

NLO_TEXT="""CT$Model = "DiracTripletEW_QCDEW_FA";
CT$GenericModel = "DiracTripletEW_QCDEW_FA";
FR$InteractionOrderPerturbativeExpansion = {QCD,QED};
R2$vertlist = {
{{F[61], F[61], A}, I*ee},
{{F[62], F[61], W}, I*ee/sw}};
UV$vertlist = {
{{F[61], F[61], A}, I*ee/FR$Eps},
{{F[63], F[62], Wbar}, I*ee/sw}};
FR$CTparam = {
MQ->MQ+FR$CT*dMQ};
"""

def evaluate(txt):
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/"ct.nlo"
        p.write_text(txt)
        return inspect(p)

class TestNLOContent(unittest.TestCase):
    def test_synthetic_one_loop_lists_structurally_complete(self):
        r=evaluate(NLO_TEXT)
        self.assertEqual(r["status"],"CONTENT_STRUCTURAL_SCREEN_PASS_NOT_PHYSICS_CERTIFIED")
        self.assertEqual(r["vertex_record_line_counts_approximate"]["UV"],2)
        self.assertEqual(r["vertex_record_line_counts_approximate"]["R2"],2)
        self.assertTrue(r["triplet_field_occurrences"]["UV"])
        self.assertEqual(sorted(r["sections_closed"]),["CTparam","R2","UV"])

    def test_missing_r2_rejected(self):
        self.assertEqual(evaluate(NLO_TEXT.replace("R2$vertlist", "wronglist"))["status"],
                         "CONTENT_REQUIRES_REVIEW")

    def test_unclosed_counterterms_rejected(self):
        self.assertEqual(evaluate(NLO_TEXT[:-3])["status"],"CONTENT_REQUIRES_REVIEW")

    def test_nonfinite_rejected(self):
        self.assertEqual(evaluate(NLO_TEXT.replace("I*ee/FR$Eps","Indeterminate"))["status"],
                         "CONTENT_REQUIRES_REVIEW")

    def test_missing_triplet_fingerprint_requires_review(self):
        s=NLO_TEXT.replace("F[61]","F[1]").replace("F[62]","F[2]").replace("F[63]","F[3]")
        self.assertEqual(evaluate(s)["status"],"CONTENT_REQUIRES_REVIEW")


if __name__=="__main__":
    unittest.main()
