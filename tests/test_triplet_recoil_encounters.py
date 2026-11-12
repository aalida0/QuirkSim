"""Tests for matched post-decay encounter differences and censored cases."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from compare_triplet_recoil_encounters import compare, quantile


def event(event_id,mode,closest_um,seed=10):
    return {
        "production_event_index":event_id,
        "params":{"mode":mode,"seed":seed,"duration_ps":2000,"lambda_ev":300},
        "initial_pdg_ids":[9000111,-9000112],
        "initial_charges_e":[1,0],
        "initial_p4_gev":[[300,10,0,100],[300,-10,0,-100]],
        "initial_charged_endpoints":1,
        "pair_pt_gev":0.,
        "n_weak_decays":1,
        "decays":[{"endpoint":0,"channel":"pion","time_ps":250}],
        "post_decay_closest_separation_cm":
            closest_um/10000 if closest_um is not None else None,
        "final_endpoint_p3_gev":[[1.,2.,3.],[4.,5.,6.]],
    }


class TestRecoilEncounters(unittest.TestCase):
    def write(self,path,rows):
        path.write_text("".join(json.dumps(row)+"\n" for row in rows))
    def test_defined_and_missing_categories(self):
        with tempfile.TemporaryDirectory() as td:
            a,b=Path(td)/"recoil.jsonl",Path(td)/"control.jsonl"
            self.write(a,[event(1,"recoil",4.),event(2,"recoil",3.),
                          event(3,"recoil",2.),event(4,"recoil",None),
                          event(5,"recoil",None)])
            self.write(b,[event(1,"no_recoil",1.),event(2,"no_recoil",5.),
                          event(3,"no_recoil",None),event(4,"no_recoil",2.),
                          event(5,"no_recoil",None)])
            out=compare(a,b)
            self.assertEqual(out["events"],5)
            self.assertEqual(out["paired_post_decay_turning_point_availability"],{
                "both_defined":2,"recoil_only_defined":1,
                "control_only_defined":1,"neither_defined":1})
            delta=out["paired_difference_recoil_minus_control_um"]
            self.assertEqual(delta["n"],2)
            self.assertAlmostEqual(delta["median"],.5)
            self.assertAlmostEqual(delta["mean"],.5)
            self.assertEqual(out["fraction_recoil_larger_conditional"],.5)
            self.assertEqual(out["fraction_recoil_smaller_conditional"],.5)
    def test_single_quantile(self):
        self.assertEqual(quantile([3.],.1),3.)
        self.assertEqual(quantile([3.],.9),3.)
    def test_reject_unmatched_seeds(self):
        with tempfile.TemporaryDirectory() as td:
            a,b=Path(td)/"a.jsonl",Path(td)/"b.jsonl"
            self.write(a,[event(1,"recoil",4.)])
            self.write(b,[event(1,"no_recoil",4.,seed=99)])
            with self.assertRaisesRegex(ValueError,"Unmatched simulation settings"):
                compare(a,b)
    def test_reject_unmatched_events(self):
        with tempfile.TemporaryDirectory() as td:
            a,b=Path(td)/"a.jsonl",Path(td)/"b.jsonl"
            self.write(a,[event(1,"recoil",4.)])
            self.write(b,[event(2,"no_recoil",4.)])
            with self.assertRaisesRegex(ValueError,"Event IDs differ"):
                compare(a,b)


if __name__=="__main__":
    unittest.main()
