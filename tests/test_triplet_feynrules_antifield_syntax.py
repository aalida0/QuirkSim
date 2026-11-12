"""Prevent silent zero-vertex FeynRules regression from wrong antifield syntax.

In FeynRules a non-self-conjugate class psi automatically defines psibar;
generic bar[psi] is not a recognized Dirac antifield. This specific bug
caused FeynmanRules[LTripletDirac] to return 0 despite three classes.
"""
from pathlib import Path
import re
import unittest

BASE=Path(__file__).resolve().parents[1]
SOURCE=BASE/"theory"/"dirac_triplet_ew"/"TripletDirac.fr"
DIAGNOSTIC=BASE/"theory"/"dirac_triplet_ew"/"00_tree_vertex_probe.wls"
STAGE1=BASE/"theory"/"dirac_triplet_ew"/"01_prepare_feynarts.wls"

class TestTripletFeynRulesSyntax(unittest.TestCase):
    def test_no_bar_function_but_automatic_feyn_rules_conjugates(self):
        source=SOURCE.read_text(encoding="utf-8")
        stripped=re.sub(r"\(\*.*?\*\)","",source,flags=re.S)
        for name in ("Qp1","Qzero","Qm2"):
            self.assertNotIn(f"bar[{name}]",stripped)
            self.assertIn(f"{name}bar",stripped)
            self.assertIn(f"ClassName -> {name}",stripped)
        self.assertIn("Mass -> {MQ, 225.0}",stripped)
        self.assertEqual(stripped.count("Mass -> MQ"),2)

    def test_runtime_guard_requires_nonzero_feynman_rules(self):
        src=DIAGNOSTIC.read_text()
        stage=STAGE1.read_text()
        self.assertIn("FR_TRIPLET_LAGRANGIAN_DECLARATION_OK",src)
        self.assertIn("FR_TREE_VERTEX_COUNT",src)
        self.assertIn("FR_STAGE1_TRIPLET_TREE_VERTEX_COUNT",stage)
        self.assertIn("Length[qsTreeVertices]<8",stage)
        self.assertIn("OnShellRenormalization[",stage)


if __name__=="__main__":
    unittest.main()
