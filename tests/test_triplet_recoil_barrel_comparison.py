"""Exact paired charge-crossing migrations on synthetic production records."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from compare_triplet_recoil_barrel import compare


def rec(event_id,mode,hit):
    rows=[{'endpoint':0,'radius_cm':3.325,
           'charge_at_crossing_e':1}] if hit else []
    return {
        'production_event_index':event_id,
        'initial_pdg_ids':[9000111,-9000112],
        'initial_charges_e':[1,0],
        'initial_p4_gev':[[250,10,0,108], [250,-10,0,-108]],
        'initial_charged_endpoints':1, 'pair_pt_gev':0,
        'params':{'lambda_ev':774,'seed':event_id,'mode':mode,
                  'b_tesla':2,'duration_ps':2000},
        'barrel_layer_radii_cm':[3.325],
        'first_outward_barrel_crossings':rows, 'n_weak_decays':1,
    }


class TestPairedMigration(unittest.TestCase):
    def write(self,path,events):
        path.write_text(''.join(json.dumps(e)+'\n' for e in events))
    def test_gain_and_loss(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=Path(tmp)/'recoil.jsonl'; b=Path(tmp)/'no_recoil.jsonl'
            self.write(a,[rec(1,'recoil',True),rec(2,'recoil',True),
                          rec(3,'recoil',False),rec(4,'recoil',False)])
            self.write(b,[rec(1,'no_recoil',False),rec(2,'no_recoil',True),
                          rec(3,'no_recoil',True),rec(4,'no_recoil',False)])
            s=compare(a,b)
            self.assertEqual(s['events'],4)
            pair=s['radii'][0]['event_level_at_least_one_charged_crossing']
            self.assertEqual(pair['both'],1)
            self.assertEqual(pair['recoil_only'],1)
            self.assertEqual(pair['no_recoil_only'],1)
            self.assertEqual(pair['neither'],1)
            self.assertEqual(pair['paired_delta_fraction'],0.0)
    def test_reject_unmatched_seed(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=Path(tmp)/'a.jsonl';b=Path(tmp)/'b.jsonl'
            self.write(a,[rec(1,'recoil',True)])
            x=rec(1,'no_recoil',True)
            x['params']['seed']=999
            self.write(b,[x])
            with self.assertRaisesRegex(ValueError,'settings other than mode'):
                compare(a,b)
    def test_reject_unmatched_sample(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=Path(tmp)/'a.jsonl';b=Path(tmp)/'b.jsonl'
            self.write(a,[rec(1,'recoil',True)])
            self.write(b,[rec(2,'no_recoil',True)])
            with self.assertRaisesRegex(ValueError,'Unmatched production'):
                compare(a,b)


if __name__=='__main__':
    unittest.main()
