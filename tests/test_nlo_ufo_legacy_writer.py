"""Regression for Python-2 parameter-card sort shipped with MG5_aMC 3.5.1."""
import ast
from pathlib import Path
import runpy
import tempfile
import unittest

patch = runpy.run_path(
    str(Path(__file__).resolve().parents[1] /
        "scripts/build_dirac_triplet_nlo_ufo.py")
)["patch_legacy_param_card_writer"]


class TestLegacyParamCardWriter(unittest.TestCase):
    def test_repairs_python2_sort_without_moving_future_import(self):
        with tempfile.TemporaryDirectory() as td:
            model = Path(td)
            writer = model / "write_param_card.py"
            writer.write_text(
                "from __future__ import absolute_import\n"
                "class Writer:\n"
                "    @staticmethod\n"
                "    def order_param(a, b):\n"
                "        return (a > b) - (a < b)\n"
                "    def ordered(self, values):\n"
                "        need_writing = list(values)\n"
                "        need_writing.sort(self.order_param)\n"
                "        return need_writing\n",
                encoding="utf-8")
            self.assertTrue(patch(model))
            content = writer.read_text(encoding="utf-8")
            ast.parse(content)
            namespace = {}
            exec(compile(content, str(writer), "exec"), namespace)
            self.assertEqual(namespace["Writer"]().ordered([3, 1, 2]), [1, 2, 3])
            self.assertFalse(patch(model))
    def test_newer_writer_untouched(self):
        with tempfile.TemporaryDirectory() as td:
            writer = Path(td) / "write_param_card.py"
            original = ("from __future__ import absolute_import\n"
                        "from functools import cmp_to_key\n"
                        "need_writing.sort(key=cmp_to_key(self.order_param))\n")
            writer.write_text(original, encoding="utf-8")
            self.assertFalse(patch(Path(td)))
            self.assertEqual(writer.read_text(encoding="utf-8"), original)


if __name__=="__main__":
    unittest.main()
