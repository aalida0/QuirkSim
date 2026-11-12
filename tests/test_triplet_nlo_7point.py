"""Verify reconstruction of 7-point NLO scale errors from saved MG5 weights."""
from pathlib import Path
import pickle
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from extract_triplet_nlo_7point import (
    get_scale_factors,read_scale_block,scale_envelope
)


class TestSevenPointScaleEnvelope(unittest.TestCase):
    def create_run(self,root):
        proc=root/"PROCESS"
        (proc/"Cards").mkdir(parents=True)
        (proc/"SubProcesses").mkdir()
        (proc/"Cards"/"run_card.dat").write_text(
            " 0.5,1.0,2.0 = rw_rscale ! MG5 reorders 1.0 first\n"
            " 2.0,0.5,1.0 = rw_fscale ! MG5 reorders 1.0 first\n")
        # Effective MG5 factors are [1,0.5,2] and [1,0.5,2] (swap 1)
        grid=[
            4.0,3.98,4.02,    # f=1, r=1,0.5,2
            4.04,4.03,4.50,   # f=0.5; r=2 excluded
            4.01,3.60,4.05,   # f=2; r=0.5 excluded
        ]
        self.assertEqual(get_scale_factors(proc/"Cards"/"run_card.dat","rw_rscale"),
                         [1.,.5,2.])
        self.assertEqual(get_scale_factors(proc/"Cards"/"run_card.dat","rw_fscale"),
                         [1.,.5,2.])
        jobs=[]
        for subnum,weight in enumerate((0.3,1.4),1):
            folder=proc/"SubProcesses"/f"P{subnum}"/f"all_G{subnum}"
            folder.mkdir(parents=True)
            jobs.append({"dirname":str(folder),"wgt_frac":weight})
            with (folder/"scale_pdf_dependence.dat").open("w") as f:
                f.write("scale variations:\n")
                f.write("-1 3 3\n")
                f.write(" ".join(f"{v/1.7:.15e}" for v in grid)+"\n")
                f.write("pdf variations:\nNNPDF40_nlo_as_01180 101\n")
                f.write(" ".join(f"{4./1.7:.15e}" for _ in range(101))+"\n")
        with (proc/"SubProcesses"/"job_status.pkl").open("wb") as f:
            pickle.dump(jobs,f)
        return proc

    def test_seven_point_excludes_antipodal_pairs(self):
        with tempfile.TemporaryDirectory() as temp:
            proc=self.create_run(Path(temp))
            res=scale_envelope(proc,expected_central=4.0)
            self.assertAlmostEqual(res["center_pb"],4.,places=11)
            self.assertAlmostEqual(res["seven_low_pb"],3.98,places=11)
            self.assertAlmostEqual(res["seven_high_pb"],4.05,places=11)
            self.assertAlmostEqual(res["seven_up_percent"],1.25,places=9)
            self.assertAlmostEqual(res["seven_down_percent"],0.5,places=9)
            self.assertAlmostEqual(res["nine_up_percent"],12.5,places=9)
            self.assertAlmostEqual(res["nine_down_percent"],10.,places=9)
            self.assertEqual(res["njobs"],2)
            self.assertEqual(len(res["points_pb"]),9)
            with self.assertRaisesRegex(ValueError,"mismatches"):
                scale_envelope(proc,expected_central=5.)

    def test_bad_scale_grid_and_missing_card(self):
        with tempfile.TemporaryDirectory() as temp:
            proc=self.create_run(Path(temp))
            file=proc/"SubProcesses"/"P1"/"all_G1"/"scale_pdf_dependence.dat"
            original=file.read_text()
            file.write_text(original.replace("-1 3 3","-1 2 2"))
            with self.assertRaisesRegex(ValueError,"Expected 9"):
                read_scale_block(file)
            file.write_text(original.replace("scale variations:","other block:"))
            with self.assertRaisesRegex(ValueError,"Missing .*scale variations"):
                read_scale_block(file)
            card=proc/"Cards"/"run_card.dat"
            card.write_text("0.5,2.0 = rw_rscale\n")
            with self.assertRaisesRegex(ValueError,"Expect exactly"):
                get_scale_factors(card,"rw_rscale")


if __name__=="__main__":
    unittest.main()
