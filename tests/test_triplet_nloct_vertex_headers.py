"""NLOCT external-leg classification, no coefficient-based false positives."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from inspect_triplet_nloct_vertex_headers import extract_first_item,inspect

SAMPLE = """CT$Model = "Triplet";
R2$vertlist = {
{{Qp1bar, Qp1, A}, coefficient_a},
{{u, ubar, Z}, function[Qp1, MQ] + coefficient_b}};
UV$vertlist = {
{{Qzerobar, Qzero, Z}, compound[300, Qp1bar, MQ]},
{{A, A}, loopfun[Qm2bar, Qp1, MQ]}};
FR$CTparam = {};
"""

def run_case(src):
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/"test.nlo"
        p.write_text(src)
        return inspect(p)

class TestNLOCTLegs(unittest.TestCase):
    def test_balanced_first_child(self):
        self.assertEqual(extract_first_item(
            b'{{Qp1bar, Qp1, A}, coeff[x, y]}'),
            b'{Qp1bar, Qp1, A}')
        self.assertEqual(extract_first_item(
            b'{{f[1,2], {a,b}, "x,y"}, 23};'),
            b'{f[1,2], {a,b}, "x,y"}')

    def test_head_isolated_from_coefficient(self):
        r=run_case(SAMPLE)
        self.assertEqual(r["status"],"HEADERS_PARSED_REQUIRE_PHYSICS_VALIDATION")
        self.assertEqual(r["records_parsed"]["UV"],2)
        self.assertEqual(r["records_parsed"]["R2"],2)
        self.assertEqual(r["explicit_triplet_external_header_counts"]["R2"],
                         {"Qp1":1,"Qp1bar":1})
        self.assertEqual(r["explicit_triplet_external_header_counts"]["UV"],
                         {"Qzero":1,"Qzerobar":1})
        self.assertNotIn("Qm2bar",r["explicit_triplet_external_header_counts"]["UV"])
        self.assertEqual(r["unparsed_examples"],[])

    def test_triplet_only_inside_coefficient_is_not_a_pass(self):
        # A huge UV expression can contain Qp1 thousands of times
        # without being a counterterm with a triplet external leg.
        sample=SAMPLE.replace("{{Qp1bar, Qp1, A}", "{{e, ebar, A}")
        sample=sample.replace("{{Qzerobar, Qzero, Z}", "{{d, dbar, Z}")
        r=run_case(sample)
        self.assertEqual(r["status"],"HEADERS_NEED_MANUAL_REVIEW")
        self.assertEqual(r["explicit_triplet_external_header_counts"]["UV"],{})
        self.assertEqual(r["explicit_triplet_external_header_counts"]["R2"],{})

    def test_broken_record_is_not_a_pass(self):
        r=run_case(SAMPLE.replace("{{u, ubar, Z}, function", "broken function"))
        self.assertEqual(r["status"],"HEADERS_NEED_MANUAL_REVIEW")
        self.assertGreater(r["records_parsed"]["R2_unparsed"],0)

if __name__=="__main__":
    unittest.main()
