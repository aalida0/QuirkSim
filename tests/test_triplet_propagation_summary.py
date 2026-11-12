"""Validate production topology and weak decay summary without simulations."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from summarize_triplet_propagation import summarize


def sample(initial, n_decay, pt):
    decays = [{'channel': 'pion','time_ps': 100+i} for i in range(n_decay)]
    return {
       'initial_charged_endpoints': initial,
       'final_charged_endpoints': initial-n_decay,
       'n_weak_decays': n_decay,
       'decays': decays,
       'pair_pt_gev': pt,
       'end_time_ps': 500.0,
    }


class TestProductionTransitions(unittest.TestCase):
    def test_user_ten_event_counts(self):
        # Recreates observed 6x(1->1), 3x(1->0), 1x(2->0)
        records = [sample(1,0,15.0)]*6 + [sample(1,1,24.0)]*3 + [sample(2,2,100.0)]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'propagation.jsonl'
            path.write_text(''.join(json.dumps(row)+'\n' for row in records))
            out = summarize(path)
            self.assertEqual(out['events'],10)
            self.assertEqual(out['initial_charged_counts'],{'1':9,'2':1})
            self.assertEqual(out['weak_decay_counts'],{'0':6,'1':3,'2':1})
            self.assertEqual(out['final_charged_counts'],{'0':4,'1':6,'2':0})
            self.assertEqual(
                {(r['initial_charged'],r['n_weak_decays'],r['final_charged']):r['events']
                 for r in out['transitions']},
                {(1,0,1):6,(1,1,0):3,(2,2,0):1})
            self.assertAlmostEqual(sum(x['fraction_all'] for x in out['transitions']),1.0)
    def test_reject_inconsistent_accounting(self):
        bad = sample(1,1,12.)
        bad['final_charged_endpoints']=1
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'bad.jsonl'
            path.write_text(json.dumps(bad)+'\n')
            with self.assertRaisesRegex(ValueError, 'mismatch'):
                summarize(path)


if __name__=='__main__':
    unittest.main()
