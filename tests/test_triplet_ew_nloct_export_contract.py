"""Check the FeynRules NLO UFO exporter receives all three NLOCT data lists.

NLOCT writes R2$vertlist, UV$vertlist, FR$CTparam.
FeynRules WriteUFO has distinct R2Vertices, UVCounterterms and CTParameters
options. Omitting CTParameters can silently discard auxiliary UV/CT coupling
objects, especially for large full-EW models.
"""
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
STAGE3=ROOT/"theory/dirac_triplet_ew"/"03_export_candidate_ufo.wls"

class TestNLOEWExportContract(unittest.TestCase):
    def test_three_nloct_lists_in_writeufo(self):
        s=STAGE3.read_text()
        self.assertIn("UVCounterterms -> UV$vertlist",s)
        self.assertIn("R2Vertices -> R2$vertlist",s)
        self.assertIn("CTParameters -> FR$CTparam",s)
        self.assertIn("!ValueQ[FR$CTparam]",s)
        self.assertIn("!ListQ[FR$CTparam]",s)
        self.assertIn("NLOCT_STAGE3_CT_PARAMETER_COUNT",s)
        self.assertIn("NLOCT_STAGE3_COUNTERTERMS_LOADED",s)

    def test_preconditions_before_full_file_get(self):
        s=STAGE3.read_text()
        for marker in (
            "stage2_log_audit.json",
            "triplet_ew_nloct_content_report.json",
            "triplet_ew_nloct_header_audit.json",
        ):
            self.assertIn(marker,s)
            self.assertLess(s.index(marker),s.index("Get[nlo]"))

    def test_candidate_deliberately_not_certified(self):
        s=STAGE3.read_text()
        self.assertIn("DiracTripletEW_QCDEW_CANDIDATE_UFO",s)
        self.assertIn("AddDecays -> False",s)
        self.assertIn("NOT CERTIFIED",s)


if __name__=="__main__":
    unittest.main()
