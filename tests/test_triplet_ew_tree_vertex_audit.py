"""Physics and negative tests for FeynRules eight-vertex FFV audit."""
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_triplet_ew_tree_vertices import audit,EXPECTED

TENSOR="Ga[Index[Lorentz, Ext[3]], Index[Spin, Ext[1]], Index[Spin, Ext[2]]]"


def make_log():
    lines=[
        "FR_TRIPLET_CLASS_DECLARATIONS_PRESENT: InputForm[{True, True, True}]",
        "FR_TRIPLET_LAGRANGIAN_DECLARATION_OK",
        "FR_TREE_PROBE_START",
        "FR_TREE_VERTEX_COUNT: 8",
    ]
    for i,(fields,(sign,factor)) in enumerate(EXPECTED.items(),1):
        stem=("I" if sign>0 else "-I")+"*"
        if factor=="ee":
            coeff=stem+"ee*"+TENSOR
        elif factor=="ee/sw":
            coeff="("+stem+"ee*"+TENSOR+")/sw"
        elif factor=="cw*ee/sw":
            coeff="("+stem+"cw*ee*"+TENSOR+")/sw"
        else:
            raise AssertionError(factor)
        fi,fo,boson=fields
        lines.append(f"FR_TREE_VERTEX_{i}: InputForm[{{{{{{{fi}, 1}}, {{{fo}, 2}}, "
                     f"{{{boson}, 3}}}}, {coeff}}}]")
    lines.append("FR_TREE_PROBE_COMPLETED_UNVALIDATED")
    return "\n".join(lines)+"\n"


def run_audit(txt):
    with tempfile.TemporaryDirectory() as d:
        path=Path(d)/"probe.log"
        path.write_text(txt)
        return audit(path)


class TestTripletVertexAudit(unittest.TestCase):
    def test_exactly_eight_expected_triplet_couplings(self):
        r=run_audit(make_log())
        self.assertTrue(r["passed"],r["missing_or_unexpected"])
        self.assertEqual(r["parsed_vertices"],8)

    def test_reject_wrong_sign_and_charged_coupling(self):
        t=make_log().replace("I*ee*"+TENSOR+")/sw",
                             "I*ee*"+TENSOR+")/(Sqrt[2]*sw)",1)
        r=run_audit(t)
        self.assertFalse(r["passed"])
        self.assertTrue(any("Unexpected coefficient" in s for s in r["missing_or_unexpected"]))

    def test_reject_missing_or_duplicate(self):
        original=make_log()
        r=run_audit(original.replace("FR_TREE_VERTEX_8:", "FR_BROKEN_VERTEX_8:",1))
        self.assertFalse(r["passed"])
        self.assertEqual(r["parsed_vertices"],7)

    def test_reject_mathematica_error(self):
        r=run_audit(make_log()+"Part::partd: unexpected Mathematica failure\n")
        self.assertFalse(r["passed"])

    def test_reject_wrong_field(self):
        r=run_audit(make_log().replace("{{Qp1bar, 1}, {Qzero, 2}, {W, 3}}",
                                      "{{Qm2bar, 1}, {Qzero, 2}, {W, 3}}"))
        self.assertFalse(r["passed"])

    def test_reject_wrong_gamma_chirality(self):
        r=run_audit(make_log().replace(TENSOR,TENSOR+"*ProjM",1))
        self.assertFalse(r["passed"])


if __name__=="__main__":
    unittest.main()
