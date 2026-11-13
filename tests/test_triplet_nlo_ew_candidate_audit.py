"""Fail-closed audit for candidate Dirac-triplet NLO QCD+EW UFOs."""
import tempfile
import json
import unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_triplet_nlo_ew_candidate import REQUIRED_FILES,inspect_ufo

class TestNLOEWGuard(unittest.TestCase):
    def make_ufo(self,path):
        path.mkdir()
        for x in REQUIRED_FILES:
            (path/x).write_text("# synthetic test fixture\n")
        # The ordinary fixture represents an exporter that wrote the optional
        # CT_parameters.py; its omission is tested explicitly below.
        (path/"CT_parameters.py").write_text("# synthetic CT parameters\n")
        (path/"particles.py").write_text("Qp1=object()\nQzero=object()\nQm2=object()\n")
        (path/"CT_vertices.py").write_text(
            "CTVertex(name='uv',type='UV',particles=[P.Qp1,P.Qp1,P.A])\n"
            "CTVertex(name='r2',type='R2',particles=[P.Qzero,P.Qp1,P.W])\n"
        )
    def test_candidate_has_triplet_uv_and_r2(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"model"
            self.make_ufo(path)
            r=inspect_ufo(path)
            self.assertTrue(r["passes_basic_triplet_ct_gate"])
            self.assertEqual(r["n_triplet_uv_vertices"],1)
            self.assertEqual(r["n_triplet_r2_vertices"],1)
            self.assertEqual(r["status"],
                             "STRUCTURAL_PREFLIGHT_ONLY_NOT_PHYSICS_VALIDATED")
    def test_no_triplet_ew_ct_is_not_nlo_ew(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"model"
            self.make_ufo(path)
            (path/"CT_vertices.py").write_text(
                "CTVertex(name='uvgluon',type='UV',particles=[P.u,P.u,P.G])\n")
            r=inspect_ufo(path)
            self.assertFalse(r["passes_basic_triplet_ct_gate"])
    def test_missing_counterterm_file_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"model"
            self.make_ufo(path)
            (path/"CT_parameters.py").unlink()
            with self.assertRaisesRegex(ValueError,"Missing mandatory"):
                inspect_ufo(path)
    def test_empty_ctparams_manifest_allows_missing_file(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"model"
            self.make_ufo(path)
            (path/"CT_parameters.py").unlink(missing_ok=True)
            (path/"quirksim_stage3_manifest.json").write_text(
                json.dumps({"schema_version":1,"ct_parameters_count":0})
            )
            self.assertTrue(inspect_ufo(path)["passes_basic_triplet_ct_gate"])

    def test_nonempty_ctparams_manifest_rejects_missing_file(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"model"
            self.make_ufo(path)
            (path/"CT_parameters.py").unlink(missing_ok=True)
            (path/"quirksim_stage3_manifest.json").write_text(
                json.dumps({"schema_version":1,"ct_parameters_count":2})
            )
            with self.assertRaisesRegex(ValueError,"CT_parameters.py"):
                inspect_ufo(path)

    def test_malformed_ctparams_manifest_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"model"
            self.make_ufo(path)
            (path/"quirksim_stage3_manifest.json").write_text(
                json.dumps({"ct_parameters_count":True})
            )
            with self.assertRaisesRegex(ValueError,"Invalid ct_parameters_count"):
                inspect_ufo(path)

    def test_missing_triplet_particles_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"model"
            self.make_ufo(path)
            (path/"particles.py").write_text("u=object()\n")
            with self.assertRaisesRegex(ValueError,"Triplet charged/neutral"):
                inspect_ufo(path)

if __name__=="__main__":
    unittest.main()
