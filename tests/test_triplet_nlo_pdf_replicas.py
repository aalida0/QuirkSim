"""Reconstruct NNPDF NLO uncertainty using MG5 job fractions and member weights."""
import math
import pickle
import statistics
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from extract_triplet_nlo_pdf_replicas import (
    read_pdf_members, pdf_replica_summary, read_central_from_log
)

class TestNNPDFReplicaExtraction(unittest.TestCase):
    def setup_jobs(self, root):
        proc = root / "process"
        (proc / "SubProcesses").mkdir(parents=True)
        jobs = []
        pdf_total = [3.989] + [3.989 + 0.01 * math.sin(i) for i in range(1, 101)]
        for job_id, fraction in ((1, 0.25), (2, 1.5)):
            folder = proc / "SubProcesses" / f"P{job_id}" / f"all_G{job_id}"
            folder.mkdir(parents=True)
            jobs.append({"dirname": str(folder), "wgt_frac": fraction})
            member_weights = [w / (0.25 + 1.5) for w in pdf_total]
            with (folder / "scale_pdf_dependence.dat").open("w") as f:
                f.write("scale variations:\n")
                f.write("10 0\n")
                f.write("1.0 1.0\n")
                f.write("pdf variations:\n")
                f.write("101 331700\n")
                f.write(" ".join(f"{v:.13e}" for v in member_weights) + "\n")
        with (proc / "SubProcesses" / "job_status.pkl").open("wb") as f:
            pickle.dump(jobs,f)
        return proc,pdf_total

    def test_replica_rms_and_central(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc, members = self.setup_jobs(Path(tmp))
            result = pdf_replica_summary(proc,expected_central=3.989)
            self.assertEqual(result["members"],101)
            self.assertEqual(result["jobs"],2)
            self.assertAlmostEqual(result["central_pb"], 3.989, places=10)
            self.assertAlmostEqual(
                result["pdf_sigma_pb"], statistics.stdev(members[1:]), places=10
            )
            with self.assertRaisesRegex(ValueError, "Wrong run"):
                pdf_replica_summary(proc,expected_central=4.5)

    def test_parser_guards_incomplete_or_corrupted_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            source=folder/"test.dat"
            source.write_text("scale variations:\n2\n1.0 2.0\n")
            with self.assertRaisesRegex(ValueError,"No 'pdf variations:'"):
                read_pdf_members(source)
            source.write_text("pdf variations:\n101\n1 2 3\n")
            with self.assertRaisesRegex(ValueError,"Expected 101"):
                read_pdf_members(source)

    def test_read_integrated_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"log"
            path.write_text(
                "Total cross section: 3.997e+00 +- 1.0e-02 pb\n"
                "Final results: Total cross section: 3.989e+00 +- 5.7e-03 pb\n"
            )
            self.assertEqual(read_central_from_log(path), 3.989)


if __name__=="__main__":
    unittest.main()
