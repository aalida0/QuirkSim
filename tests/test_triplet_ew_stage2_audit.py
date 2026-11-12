"""Check NLOCT .nlo structural audits block corruption/empty UV and R2."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_triplet_ew_stage2 import audit

LOG=("FR_STAGE1_LOG_AUDIT_PASSED_ONLY: no known rejected-vertex patterns.\n"
     "NLOCT_STAGE2_FILE_WRITTEN_UNVALIDATED: fake_output.nlo\n")
NLO=("CT$Model = \"triplet\";\nCT$GenericModel = \"triplet\";\n"
     "FR$InteractionOrderPerturbativeExpansion = {QED,QCD};\n"
     "R2$vertlist = {{F[61],F[61],A},42};\n"
     "UV$vertlist = {{F[61],F[61],A},42};\n"
     "(* this is a synthetic counterterm block, not a physical one *)\n")

def test_case(log=LOG,nlo=NLO):
    with tempfile.TemporaryDirectory() as d:
        root=Path(d)
        a,b=root/"stage2.log",root/"ct.nlo"
        a.write_text(log)
        if nlo is not None:b.write_text(nlo)
        return audit(a,b)

class TestNLOCTStage2(unittest.TestCase):
    def test_synthetic_counterterms_pass_static_preflight_only(self):
        result=test_case()
        self.assertTrue(result["preflight_pass"],result["blocking_findings"])
        self.assertEqual(result["status"],"NLOCT_STAGE2_STRUCTURAL_PREFLIGHT_PASS_NOT_PHYSICS_CERTIFIED")

    def test_empty_uv_rejected(self):
        self.assertFalse(test_case(nlo=NLO.replace(
            "UV$vertlist = {{F[61],F[61],A},42};", "UV$vertlist = {};"
        ))["preflight_pass"])

    def test_missing_r2_rejected(self):
        self.assertFalse(test_case(nlo=NLO.replace("R2$vertlist", "WRONGR2"))["preflight_pass"])

    def test_nonfinite_counterterms_rejected(self):
        self.assertFalse(test_case(nlo=NLO+"Indeterminate;")["preflight_pass"])

    def test_missing_file_rejected(self):
        self.assertFalse(test_case(nlo=None)["preflight_pass"])

    def test_math_error_rejected(self):
        self.assertFalse(test_case(log=LOG+"Part::partd: bad part\n")["preflight_pass"])

    def test_no_stage2_marker_rejected(self):
        self.assertFalse(test_case(log="FR_STAGE1_LOG_AUDIT_PASSED_ONLY\n")["preflight_pass"])

if __name__=="__main__":
    unittest.main()
