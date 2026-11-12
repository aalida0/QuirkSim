#!/usr/bin/env python3
"""Synthetic LHE tests; no MadGraph installation needed."""
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_triplet_lhe import analyze, CHANNELS


def write_lhe(file, pairs):
    lines = ['<LesHouchesEvents version="1.0">', '<init></init>']
    for a,b in pairs:
        lines.append("<event>")
        lines.append("2 1 1.0 500.0 0.007 0.1")
        for j,pid in enumerate((a,b)):
            px = 20 if j == 0 else -20
            pz = 30 if j == 0 else -30
            lines.append(f"{pid} 1 0 0 0 0 {px} 0 {pz} 230.0 225.0 0.0 9.0")
        lines.append("</event>")
    lines.append("</LesHouchesEvents>")
    content = "\n".join(lines)+"\n"
    if file.suffix == ".gz":
        with gzip.open(file, "wt") as f:
            f.write(content)
    else:
        file.write_text(content)


class TestLHEAnalyzer(unittest.TestCase):
    def test_six_channels(self):
        with tempfile.TemporaryDirectory() as td:
            file = Path(td)/"sample.lhe.gz"
            write_lhe(file, list(CHANNELS.keys()))
            result = analyze(file, sigma_pb=3.0, n_ic=2)
            self.assertEqual(result["events"],6)
            self.assertEqual([r["events"] for r in result["channels"]],[1]*6)
            self.assertEqual(result["initial_charged_endpoints"],{"0":0,"1":4,"2":2})
            self.assertEqual(result["sigma_pb_NIC_scaled"],6.0)
            self.assertAlmostEqual(sum(r["sigma_pb_NIC_scaled"] for r in result["channels"]),6.0)
            self.assertAlmostEqual(result["mean_pair_pt_GeV"],0.0)
    def test_reject_wrong_final_state(self):
        with tempfile.TemporaryDirectory() as td:
            file=Path(td)/"bad.lhe"
            write_lhe(file,[(9000112,-9000112)])
            with self.assertRaisesRegex(ValueError, "Unexpected quirk pair"):
                analyze(file)
    def test_reject_nonpositive_hidden_colors(self):
        with tempfile.TemporaryDirectory() as td:
            file=Path(td)/"small.lhe"
            write_lhe(file, [(9000111,-9000111)])
            with self.assertRaisesRegex(ValueError,"N_IC"):
                analyze(file,n_ic=0)
    def test_reject_missing_endpoints(self):
        with tempfile.TemporaryDirectory() as td:
            file=Path(td)/"bad.lhe"
            write_lhe(file,[(11,-11)])
            with self.assertRaisesRegex(ValueError,"exactly two outgoing quirks"):
                analyze(file)


if __name__=="__main__":
    unittest.main()
