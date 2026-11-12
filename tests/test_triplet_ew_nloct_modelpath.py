"""Structural guards for FeynArts' special $ModelPath resolver.

Regression: the stage-2 run failed with
  InitializeModel::noopen:
  Cannot open DiracTripletEW_QCDEW_FA/DiracTripletEW_QCDEW_FA.gen
despite the .gen file having passed the stage-1 SHA256 audit. FeynArts'
ReadModelFile uses Block[{$Path=$ModelPath}, FindFile[...]], so the
correct generated directory must be in *$ModelPath* before WriteCT.
"""
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
NLOCT=(ROOT/"theory"/"dirac_triplet_ew"/"02_calculate_nloct.wls").read_text()
DRIVER=(ROOT/"scripts"/"run_triplet_ew_stage2.sh").read_text()

class TestStage2FeynArtsModelPath(unittest.TestCase):
    def test_search_path_explicit_before_nloct(self):
        self.assertIn("$ModelPath=DeleteDuplicates[Prepend[$ModelPath,qsModelDir]]",NLOCT)
        self.assertIn("System`Private`FindFile[qsModelName<>\".gen\"]",NLOCT)
        self.assertIn("System`Private`FindFile[qsModelName<>\".mod\"]",NLOCT)
        self.assertLess(NLOCT.index("$ModelPath=DeleteDuplicates"),
                        NLOCT.index("WriteCT["))
        self.assertIn("NLOCT_FEYNARTS_INITIALIZATION_OK",NLOCT)

    def test_matching_class_and_generic_names(self):
        # Avoid erroneously passing a nested relative path through
        # FeynArts' separate $ModelPath lookup.
        self.assertIn("WriteCT[\n        qsModelName,\n        qsModelName,",NLOCT)
        self.assertNotIn('"DiracTripletEW_QCDEW_FA/DiracTripletEW_QCDEW_FA"',NLOCT)
        self.assertIn("InitializeModel[qsModelName,GenericModel->qsModelName]",NLOCT)

    def test_earlier_failed_log_preserved_and_output_not_overwritten(self):
        self.assertIn('previous_log=',DRIVER)
        self.assertIn('mv -- "$LOG" "$previous_log"',DRIVER)
        self.assertIn('if [[ -e "$NLO" ]]',DRIVER)


if __name__=="__main__":
    unittest.main()
