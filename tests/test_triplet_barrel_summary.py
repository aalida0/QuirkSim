"""Synthetic geometry-proxy unit tests, not detector efficiencies."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from summarize_triplet_barrel import summarize


def ev(charges, crossings, decay_radii=()):
    return {
        'initial_charges_e':charges,
        'barrel_layer_radii_cm':[3.325,5.05],
        'first_outward_barrel_crossings':crossings,
        'decays':[{'radius_cm':float(r),'time_ps':100.} for r in decay_radii],
    }


def crossing(i,r,charge):
    return {'endpoint':i,'radius_cm':r,'time_ps':20.,
            'z_cm':0.,'charge_at_crossing_e':charge}


class TestBarrelSummary(unittest.TestCase):
    def test_charge_neutralization_and_nonreach_distinct(self):
        # Initially charged endpoint 0 crosses IBL charged, next layer neutral.
        # Initially neutral endpoint 1 does not magically become charged.
        events=[
            ev([1,0],[crossing(0,3.325,1),crossing(0,5.05,0),
                      crossing(1,3.325,0)],(4.0,)),
            ev([1,-1],[crossing(0,3.325,1),crossing(1,3.325,-1),
                       crossing(0,5.05,1)],()),
        ]
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'r.jsonl'
            p.write_text(''.join(json.dumps(x)+'\n' for x in events))
            report=summarize(p)
        layer1,layer2=report['layers']
        self.assertEqual(report['events'],2)
        self.assertEqual(report['weak_decays_observed'],1)
        self.assertEqual(report['median_decay_radius_cm'],4.0)
        self.assertEqual(layer1['charged_crossed_while_charged'],3)
        self.assertEqual(layer1['initial_charged_endpoints'],3)
        self.assertEqual(layer1['fraction_events_two_charged_crossings'],0.5)
        self.assertEqual(layer2['charged_crossed_while_charged'],1)
        self.assertEqual(layer2['charged_crossed_after_neutralizing'],1)
        self.assertEqual(layer2['initial_charged_not_crossed_in_window'],1)
        self.assertEqual(layer2['initial_neutral_not_crossed'],1)
    def test_old_records_require_rerun(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'old.jsonl'
            p.write_text(json.dumps({'decays':[]})+'\n')
            with self.assertRaisesRegex(ValueError,'rerun'):
                summarize(p)


if __name__=='__main__':
    unittest.main()
