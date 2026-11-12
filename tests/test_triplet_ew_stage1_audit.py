"""Stage1 FeynRules/FA export must reject written but malformed .mod files."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_triplet_ew_stage1 import audit

GOOD_LOG = """FR_TRIPLET_MERGED_MODEL_LOADED
The lagrangian is hermitian.
Writing FeynArts model file into directory DiracTripletEW_QCDEW_FA
FR_STAGE1_FEYNARTS_FILES_WRITTEN_UNVALIDATED
"""

class TestStage1LogGate(unittest.TestCase):
    def check(self, suffix, models=True):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            log=root/"stage1.log"
            log.write_text(GOOD_LOG+suffix)
            m=root/"DiracTripletEW_QCDEW_FA"
            if models:
                m.mkdir()
                (m/(m.name+".mod")).write_text("fake model")
                (m/(m.name+".pars")).write_text("fake parameters")
            return audit(log,m)

    def test_clean_log_with_model_files_only_passes_log_preflight(self):
        r=self.check("")
        self.assertTrue(r["accepted_for_nloct_preflight"])
        self.assertEqual(r["status"],
                         "STAGE1_LOG_PREFLIGHT_PASS_NOT_FULL_PHYSICS_CERTIFIED")

    def test_writing_files_is_insufficient_when_mathematica_reports_errors(self):
        bad = [
            "Part::partd: Part specification Null[[1]] is longer than depth of object",
            "Set::partd: Part specification vert$[[2,1,1]] is longer",
            "Function::slotn: Slot number 1 cannot be filled",
            "StringJoin::string: String expected",
            "First::nofirst: {} has zero length",
            "Warning: The coupling V is not supported by Lorentz.gen. Vertex is ignored.",
            "Quantum number Q not conserved in vertex {W}.",
            "Warning : tadpole vertices (S), this model should be used with\n  lorentzTadpole.gen",
            "$Aborted",
        ]
        for message in bad:
            with self.subTest(message=message):
                r=self.check(message)
                self.assertFalse(r["accepted_for_nloct_preflight"])
                self.assertTrue(r["blocking_findings"])

    def test_no_model_files_must_fail(self):
        r=self.check("",models=False)
        self.assertFalse(r["accepted_for_nloct_preflight"])

    def test_no_loaded_marker_must_fail(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            log=root/"stage1.log"
            log.write_text("The lagrangian is hermitian.\n"
                           "FR_STAGE1_FEYNARTS_MODEL_OK\n")
            m=root/"DiracTripletEW_QCDEW_FA"
            m.mkdir()
            (m/(m.name+".mod")).write_text("fake")
            (m/(m.name+".pars")).write_text("fake")
            self.assertFalse(audit(log,m)["accepted_for_nloct_preflight"])


if __name__=="__main__":
    unittest.main()
